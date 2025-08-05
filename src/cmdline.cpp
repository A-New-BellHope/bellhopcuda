/*
bellhopcxx / bellhopcuda - C++/CUDA port of BELLHOP(3D) underwater acoustics simulator
Copyright (C) 2021-2023 The Regents of the University of California
Marine Physical Lab at Scripps Oceanography, c/o Jules Jaffe, jjaffe@ucsd.edu
Based on BELLHOP / BELLHOP3D, which is Copyright (C) 1983-2022 Michael B. Porter

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <https://www.gnu.org/licenses/>.
*/
#include "common_setup.hpp"

static bhc::bhcInit init;

template<bool O3D, bool R3D> int mainmain()
{
    bhc::bhcParams<O3D> params;
    bhc::bhcOutputs<O3D, R3D> outputs;
    if(!bhc::setup<O3D, R3D>(init, params, outputs)) return 1;
    if(!bhc::run<O3D, R3D>(params, outputs)) return 1;
    if(!bhc::writeout<O3D, R3D>(params, outputs, nullptr)) return 1;
    bhc::finalize<O3D, R3D>(params, outputs);
    return 0;
}

void showhelp(const char *argv0)
{
    std::cout
        << BHC_PROGRAMNAME
        " - C++/CUDA port of BELLHOP underwater acoustics simulator\n"
        "\n"
        "Copyright (C) 2021-2023 The Regents of the University of California\n"
        "Marine Physical Lab at Scripps Oceanography, c/o Jules Jaffe, jjaffe@ucsd.edu\n"
        "Based on BELLHOP, which is Copyright (C) 1983-2022 Michael B. Porter\n"
        "GPL3 licensed, no warranty, see LICENSE or https://www.gnu.org/licenses/\n"
        "\n"
        "Usage: "
        << argv0
        << " [options] FileRoot\n"
           "FileRoot is the absolute or relative path to the environment file, minus "
           "the\n"
           ".env file extension, e.g. test/in/MunkB_ray_rot .\n"
           "All command-line options may be specified with one or two dashes, e.g.\n"
           "-3 or --3 do the same thing. Furthermore, all command-line options have\n"
           "multiple synonyms which do the same thing.\n"
           "\n"
           "-?, -h, -help: Shows this help message\n"
           "-1, -singlethread: Use only one worker thread for CPU computation\n"
#if BHC_DIM_ONLY == 0
#if BHC_ENABLE_2D
           "-2, -2D: Does a 2D run. The environment file must also be 2D\n"
#endif
#if BHC_ENABLE_3D
           "-3, -3D: Does a 3D run. The environment file must also be 3D\n"
#endif
#if BHC_ENABLE_NX2D
           "-4, -Nx2D, -2D3D, -2.5D: Does a Nx2D run. The environment file must also be "
           "Nx2D\n"
#endif
#endif
           "-copy, -raycopy: Sets the behavior when there is insufficient memory to\n"
           "    allocate the requested number of full-size rays. See "
           "bhcInit::useRayCopyMode\n    in <bhc/structs.hpp> for more details\n"
#if BHC_BUILD_CUDA
           "-gpu=N, -device=N: Selects CUDA device N\n"
#endif
           "-mem=X, -memory=X: Sets the amount of memory " BHC_PROGRAMNAME
           " should use.\n"
           "    X may have a wide range of suffixes, examples: 16GiB, 8M, 100000kB\n"
           "    non-examples: 4gI, 2m, 5.3G. Default: 4GiB\n"
           "-writeenv=\"path/to/newFileRoot\": For testing purposes, writes out\n"
           "    a copy of all the input data read from the environment file etc.\n"
           "    to a new environment file and other data files. Does not run the\n"
           "    simulation\n";
}

int main(int argc, char **argv)
{
    int dimmode = BHC_DIM_ONLY;
    std::string FileRoot;
    for(int32_t i = 1; i < argc; ++i) {
        std::string s = argv[i];
        if(argv[i][0] == '-') {
            if(s.length() >= 2 && argv[i][1] == '-') { // two dashes
                s = s.substr(1);
            }
            if(s == "-1" || s == "-singlethread") {
                init.numThreads = 1;
            } else if(s == "-2" || s == "-2D") {
                dimmode = 2;
            } else if(s == "-Nx2D" || s == "-2D3D" || s == "-2.5D" || s == "-4") {
                dimmode = 4;
            } else if(s == "-3" || s == "-3D") {
                dimmode = 3;
            } else if(s == "-copy" || s == "-raycopy") {
                init.useRayCopyMode = true;
            } else if(s == "-?" || s == "-h" || s == "-help") {
                showhelp(argv[0]);
                return 0;
            } else {
                size_t equalspos = s.find("=");
                if(equalspos == std::string::npos) {
                    std::cout << "Unknown command-line option \"" << s << "\", try "
                              << argv[0] << " --help\n";
                    return 1;
                }
                std::string key   = s.substr(0, equalspos);
                std::string value = s.substr(equalspos + 1);
                if(key == "-gpu" || key == "-device") {
                    if(!bhc::isInt(value, false)) {
                        std::cout << "Value \"" << value
                                  << "\" for --gpu argument is invalid, try " << argv[0]
                                  << " --help\n";
                        return 1;
                    }
                    init.gpuIndex = std::stoi(value);
                } else if(key == "-mem" || key == "-memory") {
                    size_t multiplier = 1u;
                    size_t base       = 1000u;
                    if(bhc::endswith(value, "B") || bhc::endswith(value, "b")) {
                        value = value.substr(0, value.length() - 1);
                    }
                    if(bhc::endswith(value, "i")) {
                        base  = 1024u;
                        value = value.substr(0, value.length() - 1);
                    }
                    if(bhc::endswith(value, "k") || bhc::endswith(value, "K")) {
                        multiplier = base;
                        value      = value.substr(0, value.length() - 1);
                    } else if(bhc::endswith(value, "M")) {
                        multiplier = base * base;
                        value      = value.substr(0, value.length() - 1);
                    } else if(bhc::endswith(value, "G")) {
                        multiplier = base * base * base;
                        value      = value.substr(0, value.length() - 1);
                    }
                    if(!bhc::isInt(value, false) || (base == 1024u && multiplier == 1u)) {
                        std::cout << "Value \"" << value
                                  << "\" for --memory argument is invalid, try "
                                  << argv[0] << " --help\n";
                        return 1;
                    }
                    init.maxMemory = multiplier * std::stoi(value);
                } else {
                    std::cout << "Unknown command-line option \"-" << key << "=" << value
                              << "\", try " << argv[0] << " --help\n";
                    return 1;
                }
            }
        } else {
            if(FileRoot.empty()) {
                FileRoot = s;
            } else {
                std::cout << "Intepreting both \"" << FileRoot << "\" and \"" << s
                          << "\" as FileRoot, error\n";
                return 1;
            }
        }
    }
    if(FileRoot.empty()) {
        std::cout << "Must provide FileRoot as command-line parameter, try " << argv[0]
                  << " --help\n";
        return 1;
    }
    init.FileRoot = FileRoot.c_str();

#if BHC_DIM_ONLY > 0
    if(dimmode != BHC_DIM_ONLY) {
        std::cout << "This version of " BHC_PROGRAMNAME " was compiled to only support ";
        if constexpr(BHC_DIM_ONLY == 4) {
            std::cout << "Nx2D";
        } else {
            std::cout << BHC_DIM_ONLY << "D";
        }
        std::cout << " runs\n";
        return 1;
    }
#else
    if(dimmode < 2 || dimmode > 4) {
        std::cout << "No dimensionality specified (--2D, --Nx2D, --3D), assuming 2D\n";
        dimmode = 2;
    }
#endif

    if(dimmode == 2) {
#if BHC_ENABLE_2D
        return mainmain<false, false>();
#else
        std::cout << "This version of " BHC_PROGRAMNAME
                     " was compiled with 2D support disabled\n";
#endif
    }
    if(dimmode == 3) {
#if BHC_ENABLE_3D
        return mainmain<true, true>();
#else
        std::cout << "This version of " BHC_PROGRAMNAME
                     " was compiled with 3D support disabled\n";
#endif
    }
    if(dimmode == 4) {
#if BHC_ENABLE_NX2D
        return mainmain<true, false>();
#else
        std::cout << "This version of " BHC_PROGRAMNAME
                     " was compiled with Nx2D support disabled\n";
#endif
    }
    return 1;
}
