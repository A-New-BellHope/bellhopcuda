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

namespace bhc {

template<bool O3D> HOST_DEVICE inline int32_t GetNumJobs(
    const Position *Pos, const AnglesStructure *Angles)
{
    int32_t ret = 1;
    if(Angles->alpha.iSingle == 0) ret *= Angles->alpha.n;
    if constexpr(O3D){
        if(Angles->beta.iSingle == 0) ret *= Angles->beta.n;
        ret *= Pos->NSy;
        ret *= Pos->NSx;
    }
    ret *= Pos->NSz;
    return ret;
}

/**
 * Returns whether the job should continue.
 * `is` changed to `isrc` because `is` is used for steps
 */
template<bool O3D> HOST_DEVICE inline bool GetJobIndices(RayInitInfo &rinit, int32_t job,
    const Position *Pos, const AnglesStructure *Angles)
{
    if(Angles->alpha.iSingle >= 1){
        rinit.ialpha = Angles->alpha.iSingle - 1; // iSingle is 1-indexed because how defined in env file
    }else{
        rinit.ialpha = job % Angles->alpha.n;
        job /= Angles->alpha.n;
    }
    if constexpr(O3D){
        if(Angles->beta.iSingle >= 1){
            rinit.ibeta = Angles->beta.iSingle - 1;
        }else{
            rinit.ibeta = job % Angles->beta.n;
            job /= Angles->beta.n;
        }
        rinit.isy = job % Pos->NSy;
        job /= Pos->NSy;
        rinit.isx = job % Pos->NSx;
        job /= Pos->NSx;
    }else{
        rinit.isx = rinit.isy = rinit.ibeta = -1234567;
    }
    rinit.isz = job;
    return (rinit.isz < Pos->NSz);
}

}
