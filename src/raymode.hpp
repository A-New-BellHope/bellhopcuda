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
#include "jobs.hpp"

namespace bhc {

/**
 * Main ray tracing function for ray path output mode.
 */
template<typename CFG, bool O3D, bool R3D> HOST_DEVICE inline void MainRayMode(
    RayInitInfo &rinit, rayPt<R3D> *ray, int32_t &Nsteps, Origin<O3D, R3D> &org,
    const BdryType *ConstBdry, const BdryInfo<O3D> *bdinfo, const ReflectionInfo *refl,
    const SSPStructure *ssp, const Position *Pos, const AnglesStructure *Angles,
    const FreqInfo *freqinfo, const BeamStructure<O3D> *Beam, const BeamInfo *beaminfo,
    ErrState *errState)
{
    real DistBegTop, DistEndTop, DistBegBot, DistEndBot;
    SSPSegState iSeg;
    VEC23<O3D> xs, gradc;
    BdryState<O3D> bds;
    BdryType Bdry;

    if(!RayInit<CFG, O3D, R3D>(
           rinit, xs, ray[0], gradc, DistBegTop, DistBegBot, org, iSeg, bds, Bdry,
           ConstBdry, bdinfo, ssp, Pos, Angles, freqinfo, Beam, beaminfo, errState)) {
        Nsteps = 1;
        return;
    }

    int32_t iSmallStepCtr = 0;
    int32_t is            = 0; // index for a step along the ray

    while(true) {
        if(HasErrored(errState)) break;
        bool twoSteps = RayUpdate<CFG, O3D, R3D>(
            ray[is], ray[is + 1], ray[is + 2], DistEndTop, DistEndBot, iSmallStepCtr, org,
            iSeg, bds, Bdry, bdinfo, refl, ssp, freqinfo, Beam, xs, errState);
        if(Nsteps >= 0 && is >= Nsteps) {
            Nsteps = is + 2;
            break;
        }
        is += (twoSteps ? 2 : 1);
        if(RayTerminate<O3D, R3D>(
               ray[is], Nsteps, is, xs, iSmallStepCtr, DistBegTop, DistBegBot, DistEndTop,
               DistEndBot, org, bdinfo, Beam, errState))
            break;
    }
}

template<bool O3D, bool R3D> inline bool RunRay(
    RayInfo<O3D, R3D> *rayinfo, const bhcParams<O3D, R3D> &params, int32_t job,
    int32_t worker, RayInitInfo &rinit, int32_t &Nsteps, ErrState *errState)
{
    if(job >= rayinfo->NRays) {
        RunError(errState, BHC_ERR_JOBNUM);
        return false;
    }
    rayPt<R3D> *ray;
    if(rayinfo->isCopyMode) {
        ray = &rayinfo->WorkRayMem[worker * rayinfo->MaxPointsPerRay];
    } else {
        ray = &rayinfo->RayMem[job * rayinfo->MaxPointsPerRay];
    }
#ifdef BHC_DEBUG
    // Set to garbage values for debugging
    memset(ray, 0xFE, rayinfo->MaxPointsPerRay * sizeof(rayPt<R3D>));
#endif

    Origin<O3D, R3D> org;
    char st = params.ssp->Type;
    if(st == 'N') {
        MainRayMode<CfgSel<'R', 'G', 'N'>, O3D, R3D>(
            rinit, ray, Nsteps, org, params.Bdry, params.bdinfo, params.refl, params.ssp,
            params.Pos, params.Angles, params.freqinfo, params.Beam, params.beaminfo,
            errState);
    } else if(st == 'C') {
        MainRayMode<CfgSel<'R', 'G', 'C'>, O3D, R3D>(
            rinit, ray, Nsteps, org, params.Bdry, params.bdinfo, params.refl, params.ssp,
            params.Pos, params.Angles, params.freqinfo, params.Beam, params.beaminfo,
            errState);
    } else if(st == 'S') {
        MainRayMode<CfgSel<'R', 'G', 'S'>, O3D, R3D>(
            rinit, ray, Nsteps, org, params.Bdry, params.bdinfo, params.refl, params.ssp,
            params.Pos, params.Angles, params.freqinfo, params.Beam, params.beaminfo,
            errState);
    } else if(st == 'P') {
        MainRayMode<CfgSel<'R', 'G', 'P'>, O3D, R3D>(
            rinit, ray, Nsteps, org, params.Bdry, params.bdinfo, params.refl, params.ssp,
            params.Pos, params.Angles, params.freqinfo, params.Beam, params.beaminfo,
            errState);
    } else if(st == 'Q') {
        MainRayMode<CfgSel<'R', 'G', 'Q'>, O3D, R3D>(
            rinit, ray, Nsteps, org, params.Bdry, params.bdinfo, params.refl, params.ssp,
            params.Pos, params.Angles, params.freqinfo, params.Beam, params.beaminfo,
            errState);
    } else if(st == 'H') {
        MainRayMode<CfgSel<'R', 'G', 'H'>, O3D, R3D>(
            rinit, ray, Nsteps, org, params.Bdry, params.bdinfo, params.refl, params.ssp,
            params.Pos, params.Angles, params.freqinfo, params.Beam, params.beaminfo,
            errState);
    } else if(st == 'A') {
        MainRayMode<CfgSel<'R', 'G', 'A'>, O3D, R3D>(
            rinit, ray, Nsteps, org, params.Bdry, params.bdinfo, params.refl, params.ssp,
            params.Pos, params.Angles, params.freqinfo, params.Beam, params.beaminfo,
            errState);
    } else {
        RunError(errState, BHC_ERR_INVALID_SSP_TYPE);
        return false;
    }
    if(HasErrored(errState)) return false;

    bool ret = true;
    if(rayinfo->isCopyMode) {
        uint32_t p = AtomicFetchAdd(&rayinfo->RayMemPoints, (uint32_t)Nsteps);
        if(p + Nsteps > rayinfo->RayMemCapacity) {
            RunWarning(errState, BHC_WARN_RAYS_OUTOFMEMORY);
            rayinfo->results[job].ray = nullptr;
            ret                       = false;
        } else {
            rayinfo->results[job].ray = &rayinfo->RayMem[p];
            memcpy(rayinfo->results[job].ray, ray, Nsteps * sizeof(rayPt<R3D>));
        }
    } else {
        rayinfo->results[job].ray = ray;
    }
    rayinfo->results[job].org          = org;
    rayinfo->results[job].SrcDeclAngle = rinit.SrcDeclAngle;
    rayinfo->results[job].Nsteps       = Nsteps;

    return ret;
}

template<bool O3D, bool R3D> void InitRayMode(
    RayInfo<O3D, R3D> *rayinfo, const bhcParams<O3D, R3D> &params, uint32_t neigen);
extern template void InitRayMode<false, false>(
    RayInfo<false, false> *rayinfo, const bhcParams<false, false> &params,
    uint32_t neigen);
extern template void InitRayMode<true, false>(
    RayInfo<true, false> *rayinfo, const bhcParams<true, false> &params, uint32_t neigen);
extern template void InitRayMode<true, true>(
    RayInfo<true, true> *rayinfo, const bhcParams<true, true> &params, uint32_t neigen);

template<bool O3D, bool R3D> void PostProcessRays(
    const bhcParams<O3D, R3D> &params, RayInfo<O3D, R3D> *rayinfo);
extern template void PostProcessRays<false, false>(
    const bhcParams<false, false> &params, RayInfo<false, false> *rayinfo);
extern template void PostProcessRays<true, false>(
    const bhcParams<true, false> &params, RayInfo<true, false> *rayinfo);
extern template void PostProcessRays<true, true>(
    const bhcParams<true, true> &params, RayInfo<true, true> *rayinfo);

template<bool O3D, bool R3D> void WriteOutRays(
    const bhcParams<O3D, R3D> &params, const RayInfo<O3D, R3D> *rayinfo);
extern template void WriteOutRays<false, false>(
    const bhcParams<false, false> &params, const RayInfo<false, false> *rayinfo);
extern template void WriteOutRays<true, false>(
    const bhcParams<true, false> &params, const RayInfo<true, false> *rayinfo);
extern template void WriteOutRays<true, true>(
    const bhcParams<true, true> &params, const RayInfo<true, true> *rayinfo);

} // namespace bhc
