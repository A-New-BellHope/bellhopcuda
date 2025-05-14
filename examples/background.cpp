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

// Demonstration (and test) of non-blocking calculation.

#include <iostream>
#include <atomic>
#include <cstring>

#define BHC_DLL_IMPORT 1
#include <bhc/bhc.hpp>

void OutputCallback(const char *message)
{
    std::cout << "Out: " << message << std::endl << std::flush;
}

void PrtCallback(const char *message) { std::cout << message << std::flush; }

// Detect when the threads all complete control.
std::atomic<bool> going;
void CompletedCallback()
{
    std::cout << "CompletedCallback called.\n" << std::flush;
    going = false;
}

int main()
{
    bhc::bhcParams<true> params;
    bhc::bhcOutputs<true, true> outputs;
    bhc::bhcInit init;
    init.FileRoot          = nullptr;
    init.outputCallback    = OutputCallback;
    init.prtCallback       = PrtCallback;
    init.completedCallback = CompletedCallback;

    bhc::setup(init, params, outputs);
    // Awkward, but only rays are non-blocking.
    outputs.rayinfo->blocking = false;

    strcpy(params.Beam->RunType, "RG   3");

    // one source
    bhc::extsetup_sxsy(params, 1, 1);
    bhc::extsetup_sz(params, 1);
    params.Pos->Sx[0] = 0.0f;
    params.Pos->Sy[0] = 0.0f;
    params.Pos->Sz[0] = 1000.0f;

    // one receiver, mostly ignored in ray mode
    bhc::extsetup_rcvrbearings(params, 1);
    bhc::extsetup_rcvrranges(params, 1);
    bhc::extsetup_rcvrdepths(params, 1);
    params.Pos->RrInKm   = true;
    params.Pos->Rr[0]    = 10.0f;
    params.Pos->theta[0] = 0.0f;
    params.Pos->Rz[0]    = 1000.0f;

    // Request a bunch of rays so it takes a while.
    // May have to increase/decrease the number of rays
    //   to run in a reasonable time on your system.
    int num_rays_per_direction = 100;
    bhc::extsetup_raybearings(params, num_rays_per_direction);
    bhc::extsetup_rayelevations(params, num_rays_per_direction);
    for(int i = 0; i < num_rays_per_direction; ++i) {
        // just make the direction reasonable
        params.Angles->alpha.inDegrees = true;
        params.Angles->beta.inDegrees  = true;
        params.Angles->alpha.angles[i] = 90.0
            + 20.0 * double(i) / double(num_rays_per_direction);
        params.Angles->beta.angles[i] = 90.0
            + 20.0 * double(i) / double(num_rays_per_direction);
    }

    bhc::VEC23<true> xSSP;
    xSSP[0] = 1;
    xSSP[1] = 1;
    xSSP[2] = 10;
    float resSSP;
    std::cout << "Testing ssp call ... " << std::flush;
    if(bhc::get_ssp<true, true>(params, xSSP, resSSP)) {
        std::cout << " success " << resSSP << "\n" << std::flush;
    } else {
        std::cout << " failed " << resSSP << "\n" << std::flush;
    }

    bhc::echo(params);

    bhc::run(params, outputs);

    std::cout << "Starting the ray calculation\n" << std::flush;
    going = true;
    while(going) {
        std::cout << "   " << bhc::get_percent_progress(params) << "% done\n "
                  << std::flush;
    }
    std::cout << bhc::get_percent_progress(params) << "% threads finished\n"
              << std::flush;

    bhc::finalize(params, outputs);
    return 0;
}
