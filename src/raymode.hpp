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
#include "trace.hpp"

namespace bhc {

/**
 * Main ray tracing function for ray path output mode.
 */
HOST_DEVICE inline void MainRayMode(int32_t isrc, int32_t ialpha, real &SrcDeclAngle,
    rayPt<false> *ray, int32_t &Nsteps,
    const BdryType *ConstBdry, const BdryInfo<false> *bdinfo, const ReflectionInfo *refl,
    const SSPStructure *ssp, const Position *Pos, const AnglesStructure *Angles,
    const FreqInfo *freqinfo, const BeamStructure *Beam, const BeamInfo *beaminfo)
{
    real DistBegTop, DistEndTop, DistBegBot, DistEndBot;
    SSPSegState iSeg;
    vec2 gradc;
    BdryState<false> bds;
    BdryType Bdry;
    
    if(!RayInit(isrc, ialpha, SrcDeclAngle, ray[0], gradc, DistBegTop, DistBegBot,
        iSeg, bds, Bdry, ConstBdry, bdinfo, ssp, Pos, Angles, freqinfo, Beam, beaminfo))
    {
        Nsteps = 1;
        return;
    }
    
    int32_t iSmallStepCtr = 0;
    int32_t is = 0; // index for a step along the ray
    
    for(int32_t istep = 0; istep<MaxN-1; ++istep){
        is += RayUpdate(ray[is], ray[is+1], ray[is+2], DistEndTop, DistEndBot,
            iSmallStepCtr, iSeg, bds, Bdry, bdinfo, refl, ssp, freqinfo, Beam);
        if(RayTerminate(ray[is], Nsteps, is, DistBegTop, DistBegBot,
            DistEndTop, DistEndBot, Beam)) break;
        if(Nsteps >= 0 && is > Nsteps){
            Nsteps = is + 1;
            break;
        }
    }
}

void OpenRAYFile(LDOFile &RAYFile, std::string FileRoot, bool ThreeD, 
    const bhcParams &params);

void WriteRay2D(real alpha0, int32_t Nsteps1, LDOFile &RAYFile,
    const BdryType *Bdry, rayPt<false> *ray);
    
void InitRayMode(RayInfo *rayinfo, const bhcParams &params);
void FinalizeRayMode(RayInfo *rayinfo, std::string FileRoot, const bhcParams &params);

inline bool IsRayCopyMode(const RayInfo *rayinfo)
{
    return (size_t)rayinfo->MaxPoints < (size_t)MaxN * (size_t)rayinfo->NRays;
}

inline bool RunRay(RayInfo *rayinfo, const bhcParams &params, rayPt<false> *localmem,
    int32_t job, int32_t isrc, int32_t ialpha, int32_t &Nsteps)
{
    rayPt<false> *ray;
    if(IsRayCopyMode(rayinfo)){
        ray = localmem;
    }else{
        ray = &rayinfo->raymem[job * MaxN];
    }
    memset(ray, 0xFE, MaxN * sizeof(rayPt<false>)); //Set to garbage values for debugging
    
    real SrcDeclAngle;
    MainRayMode(isrc, ialpha, SrcDeclAngle, ray, Nsteps,
        params.Bdry, params.bdinfo, params.refl, params.ssp, params.Pos,
        params.Angles, params.freqinfo, params.Beam, params.beaminfo);
    
    bool ret = true;
    if(IsRayCopyMode(rayinfo)){
        uint32_t p = AtomicFetchAdd(&rayinfo->NPoints, (uint32_t)Nsteps);
        if(p + Nsteps > rayinfo->MaxPoints){
            std::cout << "Ran out of memory for rays\n";
            rayinfo->results[job].ray = nullptr;
            ret = false;
        }else{
            rayinfo->results[job].ray = &rayinfo->raymem[p];
            memcpy(rayinfo->results[job].ray, localmem, Nsteps * sizeof(rayPt<false>));
        }
    }else{
        rayinfo->results[job].ray = ray;
    }
    rayinfo->results[job].SrcDeclAngle = SrcDeclAngle;
    rayinfo->results[job].Nsteps = Nsteps;
    
    return ret;
}

}
