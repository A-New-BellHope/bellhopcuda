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
#pragma once
#include "common.hpp"
#include "raymode.hpp"
#include "tlmode.hpp"
#include "eigenrays.hpp"
#include "arrivals.hpp"

namespace bhc {

HOST_DEVICE inline int32_t GetNumJobs(const Position *Pos, const AnglesStructure *Angles)
{
    return Pos->NSz * (Angles->iSingle_alpha >= 1 ? 1 : Angles->Nalpha);
}

/**
 * Returns whether the job should continue.
 * `is` changed to `isrc` because `is` is used for steps
 */
HOST_DEVICE inline bool GetJobIndices(RayInitInfo &rinit, int32_t job,
    const Position *Pos, const AnglesStructure *Angles)
{
    if(Angles->iSingle_alpha >= 1){
        rinit.isz = job;
        rinit.ialpha = Angles->iSingle_alpha - 1; // iSingle_alpha is 1-indexed because how defined in env file
    }else{
        rinit.isz = job / Angles->Nalpha;
        rinit.ialpha = job % Angles->Nalpha;
    }
    rinit.isx = rinit.isy = rinit.ibeta = -1234567; // TODO
    return (rinit.isz < Pos->NSz);
}

inline void InitSelectedMode(const bhcParams &params, bhcOutputs &outputs, bool singlethread)
{
    // Common
    // irregular or rectilinear grid
    params.Pos->NRz_per_range = (params.Beam->RunType[4] == 'I') ? 1 : params.Pos->NRz;
    
    // Mode specific
    if(params.Beam->RunType[0] == 'R'){
        InitRayMode(outputs.rayinfo, params);
    }else if(params.Beam->RunType[0] == 'C' || params.Beam->RunType[0] == 'S' || params.Beam->RunType[0] == 'I'){
        InitTLMode(outputs.uAllSources, params.Pos);
    }else if(params.Beam->RunType[0] == 'E'){
        InitEigenMode(outputs.eigen);
    }else if(params.Beam->RunType[0] == 'A' || params.Beam->RunType[0] == 'a'){
        InitArrivalsMode(outputs.arrinfo, singlethread, params.Pos, *(PrintFileEmu*)params.internal);
    }else{
        std::cout << "Invalid RunType " << params.Beam->RunType[0] << "\n";
        std::abort();
    }
}

bool run_cxx(const bhcParams &params, bhcOutputs &outputs, bool singlethread);

}
