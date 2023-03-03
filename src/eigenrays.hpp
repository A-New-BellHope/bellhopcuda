/*
bellhopcxx / bellhopcuda - C++/CUDA port of BELLHOP underwater acoustics simulator
Copyright (C) 2021-2023 The Regents of the University of California
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

namespace bhc {

HOST_DEVICE inline void RecordEigenHit(
    int32_t itheta, int32_t ir, int32_t iz, int32_t is, const RayInitInfo &rinit,
    EigenInfo *eigen)
{
    uint32_t mi = AtomicFetchAdd(&eigen->neigen, 1u);
    if(mi >= eigen->memsize) return;
    // printf("Eigenray hit %d ir %d iz %d isrc %d ialpha %d is %d\n",
    //     mi, ir, iz, isrc, ialpha, is);
    eigen->hits[mi].is     = is;
    eigen->hits[mi].iz     = iz;
    eigen->hits[mi].ir     = ir;
    eigen->hits[mi].itheta = itheta;
    eigen->hits[mi].isx    = rinit.isx;
    eigen->hits[mi].isy    = rinit.isy;
    eigen->hits[mi].isz    = rinit.isz;
    eigen->hits[mi].ialpha = rinit.ialpha;
    eigen->hits[mi].ibeta  = rinit.ibeta;
}

template<bool O3D, bool R3D> inline void InitEigenMode(
    EigenInfo *eigen, const bhcParams<O3D, R3D> &params)
{
#warning TODO InitEigenMode
    constexpr uint32_t maxhits = 1000000u;
    eigen->neigen              = 0;
    eigen->memsize             = maxhits;
    trackallocate(params, "eigenray hits", eigen->hits, maxhits);
}

template<bool O3D, bool R3D> void PostProcessEigenrays(
    const bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs);
extern template void PostProcessEigenrays<false, false>(
    const bhcParams<false, false> &params, bhcOutputs<false, false> &outputs);
extern template void PostProcessEigenrays<true, false>(
    const bhcParams<true, false> &params, bhcOutputs<true, false> &outputs);
extern template void PostProcessEigenrays<true, true>(
    const bhcParams<true, true> &params, bhcOutputs<true, true> &outputs);

template<bool O3D, bool R3D> inline void WriteOutEigenrays(
    const bhcParams<O3D, R3D> &params, const bhcOutputs<O3D, R3D> &outputs)
{
    WriteOutRays<O3D, R3D>(params, outputs.rayinfo);
}

} // namespace bhc
