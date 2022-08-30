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

template<bool O3D, bool R3D> inline void InitSelectedMode(
    const bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs, bool singlethread)
{
    // Common
    int32_t ns = params.Pos->NSx * params.Pos->NSy * params.Pos->NSz;
    BASSERT(ns >= 1);
    
    // irregular or rectilinear grid
    params.Pos->NRz_per_range = (params.Beam->RunType[4] == 'I') ? 1 : params.Pos->NRz;
    
    // Mode specific
    if(params.Beam->RunType[0] == 'R'){
        InitRayMode<O3D, R3D>(outputs.rayinfo, params);
    }else if(IsTLRun(params.Beam)){
        InitTLMode(outputs.uAllSources, params.Pos);
    }else if(params.Beam->RunType[0] == 'E'){
        InitEigenMode(outputs.eigen);
    }else if(IsArrivalsRun(params.Beam)){
        InitArrivalsMode(outputs.arrinfo, singlethread, R3D, params.Pos,
            *(PrintFileEmu*)params.internal);
    }else{
        std::cout << "Invalid RunType " << params.Beam->RunType[0] << "\n";
        std::abort();
    }
}

template<bool O3D, bool R3D> bool run_cxx(
    const bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs, bool singlethread);
extern template bool run_cxx<false, false>(
    const bhcParams<false, false> &params, bhcOutputs<false, false> &outputs, bool singlethread);
extern template bool run_cxx<true, false>(
    const bhcParams<true, false> &params, bhcOutputs<true, false> &outputs, bool singlethread);
extern template bool run_cxx<true, true>(
    const bhcParams<true, true> &params, bhcOutputs<true, true> &outputs, bool singlethread);


}
