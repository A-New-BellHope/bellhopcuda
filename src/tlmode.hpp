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

template<bool O3D, bool R3D> void PostProcessTL(
    const bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs);
extern template void PostProcessTL<false, false>(
    const bhcParams<false, false> &params, bhcOutputs<false, false> &outputs);
extern template void PostProcessTL<true, false>(
    const bhcParams<true, false> &params, bhcOutputs<true, false> &outputs);
extern template void PostProcessTL<true, true>(
    const bhcParams<true, true> &params, bhcOutputs<true, true> &outputs);

template<bool O3D, bool R3D> void WriteOutTL(
    const bhcParams<O3D, R3D> &params, const bhcOutputs<O3D, R3D> &outputs);
extern template void WriteOutTL<false, false>(
    const bhcParams<false, false> &params, const bhcOutputs<false, false> &outputs);
extern template void WriteOutTL<true, false>(
    const bhcParams<true, false> &params, const bhcOutputs<true, false> &outputs);
extern template void WriteOutTL<true, true>(
    const bhcParams<true, true> &params, const bhcOutputs<true, true> &outputs);

/**
 * Main ray tracing function for TL, eigen, and arrivals runs.
 */
template<typename CFG, bool O3D, bool R3D> HOST_DEVICE inline void MainFieldModes(
    RayInitInfo &rinit, cpxf *uAllSources, const BdryType *ConstBdry,
    const BdryInfo<O3D> *bdinfo, const ReflectionInfo *refl, const SSPStructure *ssp,
    const Position *Pos, const AnglesStructure *Angles, const FreqInfo *freqinfo,
    const BeamStructure<O3D> *Beam, const BeamInfo *beaminfo, EigenInfo *eigen,
    const ArrInfo *arrinfo, ErrState *errState)
{
    real DistBegTop, DistEndTop, DistBegBot, DistEndBot;
    SSPSegState iSeg;
    VEC23<O3D> xs, gradc;
    BdryState<O3D> bds;
    BdryType Bdry;
    Origin<O3D, R3D> org;

    rayPt<R3D> point0, point1, point2;
    point2.c = NAN; // Silence incorrect g++ warning about maybe uninitialized;
    // it is always set when doing two steps, and not used otherwise
    InfluenceRayInfo<R3D> inflray;

    if(!RayInit<CFG, O3D, R3D>(
           rinit, xs, point0, gradc, DistBegTop, DistBegBot, org, iSeg, bds, Bdry,
           ConstBdry, bdinfo, ssp, Pos, Angles, freqinfo, Beam, beaminfo, errState)) {
        return;
    }

    Init_Influence<CFG, O3D, R3D>(
        inflray, point0, rinit, gradc, Pos, org, ssp, iSeg, Angles, freqinfo, Beam,
        errState);

    int32_t iSmallStepCtr = 0;
    int32_t is            = 0; // index for a step along the ray
    int32_t Nsteps        = 0; // not actually needed in TL mode, debugging only

    for(int32_t istep = 0; istep < MaxN - 1; ++istep) {
        if(HasErrored(errState)) break;
        bool twoSteps = RayUpdate<CFG, O3D, R3D>(
            point0, point1, point2, DistEndTop, DistEndBot, iSmallStepCtr, org, iSeg, bds,
            Bdry, bdinfo, refl, ssp, freqinfo, Beam, xs, errState);
        if(!Step_Influence<CFG, O3D, R3D>(
               point0, point1, inflray, is, uAllSources, ConstBdry, org, ssp, iSeg, Pos,
               Beam, eigen, arrinfo, errState)) {
#ifdef STEP_DEBUGGING
            printf("Step_Influence terminated ray\n");
#endif
            break;
        }
        ++is;
        if(twoSteps) {
            if(!Step_Influence<CFG, O3D, R3D>(
                   point1, point2, inflray, is, uAllSources, ConstBdry, org, ssp, iSeg,
                   Pos, Beam, eigen, arrinfo, errState))
                break;
            point0 = point2;
            ++is;
        } else {
            point0 = point1;
        }
        if(RayTerminate<O3D, R3D>(
               point0, Nsteps, is, xs, iSmallStepCtr, DistBegTop, DistBegBot, DistEndTop,
               DistEndBot, org, bdinfo, Beam, errState))
            break;
    }

    // printf("Nsteps %d\n", Nsteps);
}

} // namespace bhc
