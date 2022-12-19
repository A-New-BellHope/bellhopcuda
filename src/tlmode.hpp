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
    size_t n = (size_t)Pos->NSz * (size_t)Pos->NSx * (size_t)Pos->NSy
        * (size_t)Pos->Ntheta * (size_t)Pos->NRz_per_range * (size_t)Pos->NRr;
    checkallocate(uAllSources, n);
    memset(uAllSources, 0, n * sizeof(cpxf));
}

template<bool O3D, bool R3D> void FinalizeTLMode(
    std::string FileRoot, const bhcParams<O3D, R3D> &params,
    bhcOutputs<O3D, R3D> &outputs);
extern template void FinalizeTLMode<false, false>(
    std::string FileRoot, const bhcParams<false, false> &params,
    bhcOutputs<false, false> &outputs);
extern template void FinalizeTLMode<true, false>(
    std::string FileRoot, const bhcParams<true, false> &params,
    bhcOutputs<true, false> &outputs);
extern template void FinalizeTLMode<true, true>(
    std::string FileRoot, const bhcParams<true, true> &params,
    bhcOutputs<true, true> &outputs);

/**
 * Main ray tracing function for TL, eigen, and arrivals runs.
 */
template<bool O3D, bool R3D> HOST_DEVICE inline void MainFieldModes(
    RayInitInfo &rinit, cpxf *uAllSources, const BdryType *ConstBdry,
    const BdryInfo<O3D> *bdinfo, const ReflectionInfo *refl, const SSPStructure *ssp,
    const Position *Pos, const AnglesStructure *Angles, const FreqInfo *freqinfo,
    const BeamStructure<O3D> *Beam, const BeamInfo *beaminfo, EigenInfo *eigen,
    const ArrInfo *arrinfo)
{
    real DistBegTop, DistEndTop, DistBegBot, DistEndBot;
    SSPSegState iSeg;
    VEC23<O3D> xs, gradc;
    BdryState<O3D> bds;
    BdryType Bdry;
    Origin<O3D, R3D> org;

    rayPt<R3D> point0, point1, point2;
    InfluenceRayInfo<R3D> inflray;

    if(!RayInit<O3D, R3D>(
           rinit, xs, point0, gradc, DistBegTop, DistBegBot, org, iSeg, bds, Bdry,
           ConstBdry, bdinfo, ssp, Pos, Angles, freqinfo, Beam, beaminfo)) {
        return;
    }

    Init_Influence<O3D, R3D>(
        inflray, point0, rinit, gradc, Pos, org, ssp, iSeg, Angles, freqinfo, Beam);

    int32_t iSmallStepCtr = 0;
    int32_t is            = 0; // index for a step along the ray
    int32_t Nsteps        = 0; // not actually needed in TL mode, debugging only

    for(int32_t istep = 0; istep < MaxN - 1; ++istep) {
        int32_t dStep = RayUpdate<O3D, R3D>(
            point0, point1, point2, DistEndTop, DistEndBot, iSmallStepCtr, org, iSeg, bds,
            Bdry, bdinfo, refl, ssp, freqinfo, Beam, xs);
        if(!Step_Influence<O3D, R3D>(
               point0, point1, inflray, is, uAllSources, ConstBdry, org, ssp, iSeg, Pos,
               Beam, eigen, arrinfo)) {
#ifdef STEP_DEBUGGING
            GlobalLog("Step_Influence terminated ray\n");
#endif
            break;
        }
        ++is;
        if(dStep == 2) {
            if(!Step_Influence<O3D, R3D>(
                   point1, point2, inflray, is, uAllSources, ConstBdry, org, ssp, iSeg,
                   Pos, Beam, eigen, arrinfo))
                break;
            point0 = point2;
            ++is;
        } else if(dStep == 1) {
            point0 = point1;
        } else {
            GlobalLog("Invalid dStep: %d\n", dStep);
            bail();
        }
        if(RayTerminate<O3D, R3D>(
               point0, Nsteps, is, xs, iSmallStepCtr, DistBegTop, DistBegBot, DistEndTop,
               DistEndBot, org, bdinfo, Beam))
            break;
    }

    // GlobalLog("Nsteps %d\n", Nsteps);
}

} // namespace bhc
