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
#pragma once
#include "common_run.hpp"

namespace bhc {

HOST_DEVICE inline void RecordEigenHit(
    int32_t itheta, int32_t ir, int32_t iz, int32_t is, const RayInitInfo &rinit,
    EigenInfo *eigen)
{
    int32_t mi = AtomicFetchAdd(&eigen->neigen, 1);
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

} // namespace bhc
