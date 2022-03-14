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

HOST_DEVICE inline void RecordEigenHit(int32_t ir, int32_t iz, 
    int32_t isrc, int32_t ialpha, int32_t is, EigenInfo *eigen)
{
    uint32_t mi = AtomicFetchAdd(&eigen->neigen, 1u);
    if(mi >= eigen->memsize) return;
    // printf("Eigenray hit %d ir %d iz %d isrc %d ialpha %d is %d\n",
    //     mi, ir, iz, isrc, ialpha, is);
    eigen->hits[mi].ir = ir;
    eigen->hits[mi].iz = iz;
    eigen->hits[mi].isrc = isrc;
    eigen->hits[mi].ialpha = ialpha;
    eigen->hits[mi].is = is;
}

inline void InitEigenMode(EigenInfo *eigen)
{
    constexpr uint32_t maxhits = 1000000u;
    eigen->neigen = 0;
    eigen->memsize = maxhits;
    checkallocate(eigen->hits, maxhits);
}

void FinalizeEigenMode(const bhcParams &params, bhcOutputs &outputs, 
    std::string FileRoot, bool singlethread);

}
