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
#include "ssp.hpp"
#include "influence.hpp"

namespace bhc {

/**
 * for a TL calculation, allocate space for the pressure matrix
 */
inline void InitTLMode(cpxf *&uAllSources, const Position *Pos)
{
    size_t n = Pos->NSz * Pos->NRz_per_range * Pos->NRr;
    checkallocate(uAllSources, n);
    memset(uAllSources, 0, n * sizeof(cpxf));
}

template<bool O3D, bool R3D> void FinalizeTLMode(
    std::string FileRoot, const bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs);
extern template void FinalizeTLMode<false, false>(
    std::string FileRoot, const bhcParams<false, false> &params, bhcOutputs<false, false> &outputs);
extern template void FinalizeTLMode<true, false>(
    std::string FileRoot, const bhcParams<true, false> &params, bhcOutputs<true, false> &outputs);
extern template void FinalizeTLMode<true, true>(
    std::string FileRoot, const bhcParams<true, true> &params, bhcOutputs<true, true> &outputs);

/**
 * Main ray tracing function for TL, eigen, and arrivals runs.
 */
HOST_DEVICE inline void MainFieldModes(RayInitInfo &rinit, cpxf *uAllSources,
    const BdryType *ConstBdry, const BdryInfo<false> *bdinfo, const ReflectionInfo *refl,
    const SSPStructure *ssp, const Position *Pos, const AnglesStructure *Angles,
    const FreqInfo *freqinfo, const BeamStructure *Beam, const BeamInfo *beaminfo,
    EigenInfo *eigen, const ArrInfo *arrinfo)
{
    real DistBegTop, DistEndTop, DistBegBot, DistEndBot;
    SSPSegState iSeg;
    vec2 xs, gradc;
    BdryState<false> bds;
    BdryType Bdry;
    Origin<false, false> org;
    
    rayPt<false> point0, point1, point2;
    InfluenceRayInfo inflray;
    
    if(!RayInit<false, false>(rinit, xs, point0, gradc, DistBegTop, DistBegBot,
        org, iSeg, bds, Bdry, ConstBdry, bdinfo, ssp, Pos, Angles, freqinfo, Beam, beaminfo))
    {
        return;
    }
    
    Init_Influence(inflray, point0, rinit, gradc,
        Pos, org, ssp, iSeg, Angles, freqinfo, Beam);
    
    cpxf *u = &uAllSources[rinit.isz * Pos->NRz_per_range * Pos->NRr];
    int32_t iSmallStepCtr = 0;
    int32_t is = 0; // index for a step along the ray
    int32_t Nsteps = 0; // not actually needed in TL mode, debugging only
    
    for(int32_t istep = 0; istep<MaxN-1; ++istep){
        int32_t dStep = RayUpdate<false, false>(point0, point1, point2,
            DistEndTop, DistEndBot, iSmallStepCtr,
            org, iSeg, bds, Bdry, bdinfo, refl, ssp, freqinfo, Beam);
        if(!Step_Influence(point0, point1, inflray, is, u, 
            ConstBdry, org, ssp, iSeg, Pos, Beam, eigen, arrinfo)){
            //printf("Step_Influence terminated ray\n");
            break;
        }
        ++is;
        if(dStep == 2){
            if(!Step_Influence(point1, point2, inflray, is, u, 
                ConstBdry, org, ssp, iSeg, Pos, Beam, eigen, arrinfo)) break;
            point0 = point2;
            ++is;
        }else if(dStep == 1){
            point0 = point1;
        }else{
            printf("Invalid dStep: %d\n", dStep);
            bail();
        }
        if(RayTerminate(point0, Nsteps, is, xs, iSmallStepCtr,
            DistBegTop, DistBegBot, DistEndTop, DistEndBot, org, bdinfo, Beam)) break;
    }
    
    //printf("Nsteps %d\n", Nsteps);
}

}
