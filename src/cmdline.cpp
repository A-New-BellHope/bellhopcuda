/*
bellhopcxx / bellhopcuda - C++/CUDA port of BELLHOP underwater acoustics simulator
Copyright (C) 2021-2022 The Regents of the University of California
c/o Jules Jaffe team at SIO / UCSD, jjaffe@ucsd.edu
Based on BELLHOP, which is Copyright (C) 1983-2020 Michael B. Porter

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
#include "common.hpp"
#include "raymode.hpp"
#include "tlmode.hpp"
#include "eigenrays.hpp"
#include "arrivals.hpp"

static std::string FileRoot;
static bool singlethread = false;

template<bool O3D, bool R3D> int mainmain()
{
    bhc::bhcParams<O3D, R3D> params;
    bhc::bhcOutputs<O3D, R3D> outputs;
    if(!bhc::setup<O3D, R3D>(FileRoot.c_str(), nullptr, params, outputs)) return 1;

    bhc::Stopwatch sw;
    sw.tick();
    if(!bhc::run<O3D, R3D>(params, outputs, singlethread)) return 1;
    sw.tock();

    if(IsRayRun(params.Beam)) {
        // Ray mode
        bhc::FinalizeRayMode<O3D, R3D>(outputs.rayinfo, FileRoot, params);
    } else if(IsTLRun(params.Beam)) {
        // TL mode
        bhc::FinalizeTLMode(FileRoot, params, outputs);
    } else if(IsEigenraysRun(params.Beam)) {
        // Eigenrays mode
        bhc::FinalizeEigenMode<O3D, R3D>(params, outputs, FileRoot, singlethread);
    } else if(IsArrivalsRun(params.Beam)) {
        // Arrivals mode
        bhc::FinalizeArrivalsMode<O3D, R3D>(
            outputs.arrinfo, params.Pos, params.freqinfo, params.Beam, FileRoot);
    } else {
        std::cout << "Invalid RunType " << params.Beam->RunType[0] << "\n";
        std::abort();
    }

    bhc::finalize<O3D, R3D>(params, outputs);
    return 0;
}

int main(int argc, char **argv)
{
    int dimmode = BHC_DIMMODE;
    for(int32_t i = 1; i < argc; ++i) {
        std::string s = argv[i];
        if(argv[i][0] == '-') {
            if(s.length() >= 2 && argv[i][1] == '-') { // two dashes
                s = s.substr(1);
            }
            if(s == "-1" || s == "-singlethread") {
                singlethread = true;
            } else if(s == "-2" || s == "-2D") {
                dimmode = 2;
            } else if(s == "-Nx2D" || s == "-2D3D" || s == "-2.5D" || s == "-4") {
                dimmode = 4;
            } else if(s == "-3" || s == "-3D") {
                dimmode = 3;
            } else {
                std::cout << "Unknown command-line option \"" << s << "\"\n";
                std::abort();
            }
        } else {
            if(FileRoot.empty()) {
                FileRoot = s;
            } else {
                std::cout << "Intepreting both \"" << FileRoot << "\" and \"" << s
                          << "\" as FileRoot, error\n";
                std::abort();
            }
        }
    }
    if(FileRoot.empty()) {
        std::cout << "Must provide FileRoot as command-line parameter\n";
        std::abort();
    }
#if BHC_DIMMODE
    if(dimmode != BHC_DIMMODE) {
        std::cout << "Cannot change dimensionality, this is " BHC_PROGRAMNAME "\n";
        std::abort();
    }
#else
    if(dimmode == 0) {
        std::cout << "No dimensionality specified (--2D, --Nx2D, --3D), assuming 2D\n";
        dimmode = 2;
    }
#endif

#if BHC_ENABLE_2D
    if(dimmode == 2) { return mainmain<false, false>(); }
#endif
#if BHC_ENABLE_3D
    if(dimmode == 3) { return mainmain<true, true>(); }
#endif
#if BHC_ENABLE_NX2D
    if(dimmode == 4) { return mainmain<true, false>(); }
#endif
    std::cout << "Internal error\n";
    std::abort();
}
