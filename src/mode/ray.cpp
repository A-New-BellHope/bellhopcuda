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
#include "ray.hpp"
#include "../trace.hpp"

#include <vector>

namespace bhc { namespace mode {

template<bool O3D, bool R3D> bool RunRay(
    RayInfo<O3D, R3D> *rayinfo, const bhcParams<O3D> &params, int32_t job, int32_t worker,
    RayInitInfo &rinit, int32_t &Nsteps, ErrState *errState)
{
    if(job >= rayinfo->NRays || worker >= GetInternal(params)->numThreads) {
        RunError(errState, BHC_ERR_JOBNUM);
        return false;
    }
    rayPt<R3D> *ray;
    if(rayinfo->isCopyMode) {
        ray = &rayinfo->WorkRayMem[worker * rayinfo->MaxPointsPerRay];
    } else {
        ray = &rayinfo->RayMem[(size_t)job * rayinfo->MaxPointsPerRay];
    }
#ifdef BHC_DEBUG
    // Set to garbage values for debugging
    memset(ray, 0xFE, rayinfo->MaxPointsPerRay * sizeof(rayPt<R3D>));
#endif

    Origin<O3D, R3D> org;
    char st = params.ssp->Type;
    if(st == 'N') {
        MainRayMode<CfgSel<'R', 'G', 'N'>, O3D, R3D>(
            rinit, ray, Nsteps, rayinfo->MaxPointsPerRay, org, params.Bdry, params.bdinfo,
            params.refl, params.ssp, params.Pos, params.Angles, params.freqinfo,
            params.Beam, params.sbp, errState);
    } else if(st == 'C') {
        MainRayMode<CfgSel<'R', 'G', 'C'>, O3D, R3D>(
            rinit, ray, Nsteps, rayinfo->MaxPointsPerRay, org, params.Bdry, params.bdinfo,
            params.refl, params.ssp, params.Pos, params.Angles, params.freqinfo,
            params.Beam, params.sbp, errState);
    } else if(st == 'S') {
        MainRayMode<CfgSel<'R', 'G', 'S'>, O3D, R3D>(
            rinit, ray, Nsteps, rayinfo->MaxPointsPerRay, org, params.Bdry, params.bdinfo,
            params.refl, params.ssp, params.Pos, params.Angles, params.freqinfo,
            params.Beam, params.sbp, errState);
    } else if(st == 'P') {
        MainRayMode<CfgSel<'R', 'G', 'P'>, O3D, R3D>(
            rinit, ray, Nsteps, rayinfo->MaxPointsPerRay, org, params.Bdry, params.bdinfo,
            params.refl, params.ssp, params.Pos, params.Angles, params.freqinfo,
            params.Beam, params.sbp, errState);
    } else if(st == 'Q') {
        MainRayMode<CfgSel<'R', 'G', 'Q'>, O3D, R3D>(
            rinit, ray, Nsteps, rayinfo->MaxPointsPerRay, org, params.Bdry, params.bdinfo,
            params.refl, params.ssp, params.Pos, params.Angles, params.freqinfo,
            params.Beam, params.sbp, errState);
    } else if(st == 'H') {
        MainRayMode<CfgSel<'R', 'G', 'H'>, O3D, R3D>(
            rinit, ray, Nsteps, rayinfo->MaxPointsPerRay, org, params.Bdry, params.bdinfo,
            params.refl, params.ssp, params.Pos, params.Angles, params.freqinfo,
            params.Beam, params.sbp, errState);
    } else if(st == 'A') {
        MainRayMode<CfgSel<'R', 'G', 'A'>, O3D, R3D>(
            rinit, ray, Nsteps, rayinfo->MaxPointsPerRay, org, params.Bdry, params.bdinfo,
            params.refl, params.ssp, params.Pos, params.Angles, params.freqinfo,
            params.Beam, params.sbp, errState);
    } else {
        RunError(errState, BHC_ERR_INVALID_SSP_TYPE);
        return false;
    }
    if(HasErrored(errState)) return false;

    bool ret = true;
    if(rayinfo->isCopyMode) {
        size_t p = AtomicFetchAdd(&rayinfo->RayMemPoints, (size_t)Nsteps);
        if(p + (size_t)Nsteps > rayinfo->RayMemCapacity) {
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

#if BHC_ENABLE_2D
template bool RunRay<false, false>(
    RayInfo<false, false> *rayinfo, const bhcParams<false> &params, int32_t job,
    int32_t worker, RayInitInfo &rinit, int32_t &Nsteps, ErrState *errState);
#endif
#if BHC_ENABLE_NX2D
template bool RunRay<true, false>(
    RayInfo<true, false> *rayinfo, const bhcParams<true> &params, int32_t job,
    int32_t worker, RayInitInfo &rinit, int32_t &Nsteps, ErrState *errState);
#endif
#if BHC_ENABLE_3D
template bool RunRay<true, true>(
    RayInfo<true, true> *rayinfo, const bhcParams<true> &params, int32_t job,
    int32_t worker, RayInitInfo &rinit, int32_t &Nsteps, ErrState *errState);
#endif

template<bool O3D, bool R3D> void RayModeWorker(
    const bhcParams<O3D> &params, bhcOutputs<O3D, R3D> &outputs, int32_t worker,
    ErrState *errState)
{
    SetupThread();
    while(true) {
        int32_t job    = GetInternal(params)->sharedJobID++;
        int32_t Nsteps = -1;
        RayInitInfo rinit;
        if(!GetJobIndices<O3D>(rinit, job, params.Pos, params.Angles)) break;
        if(!RunRay<O3D, R3D>(
               outputs.rayinfo, params, job, worker, rinit, Nsteps, errState)) {
            break;
        }
    }
}

#if BHC_ENABLE_2D
template void RayModeWorker<false, false>(
    const bhcParams<false> &params, bhcOutputs<false, false> &outputs, int32_t worker,
    ErrState *errState);
#endif
#if BHC_ENABLE_NX2D
template void RayModeWorker<true, false>(
    const bhcParams<true> &params, bhcOutputs<true, false> &outputs, int32_t worker,
    ErrState *errState);
#endif
#if BHC_ENABLE_3D
template void RayModeWorker<true, true>(
    const bhcParams<true> &params, bhcOutputs<true, true> &outputs, int32_t worker,
    ErrState *errState);
#endif

template<bool O3D, bool R3D> void RunRayMode(
    bhcParams<O3D> &params, bhcOutputs<O3D, R3D> &outputs)
{
    ErrState errState;
    ResetErrState(&errState);
    GetInternal(params)->sharedJobID = 0;
    int32_t numThreads               = GetInternal(params)->numThreads;
    std::vector<std::thread> threads;
    for(int32_t i = 0; i < numThreads; ++i)
        threads.push_back(std::thread(
            RayModeWorker<O3D, R3D>, std::ref(params), std::ref(outputs), i, &errState));
    for(int32_t i = 0; i < numThreads; ++i) threads[i].join();
    CheckReportErrors(GetInternal(params), &errState);
}

#if BHC_ENABLE_2D
template void RunRayMode<false, false>(
    bhcParams<false> &params, bhcOutputs<false, false> &outputs);
#endif
#if BHC_ENABLE_NX2D
template void RunRayMode<true, false>(
    bhcParams<true> &params, bhcOutputs<true, false> &outputs);
#endif
#if BHC_ENABLE_3D
template void RunRayMode<true, true>(
    bhcParams<true> &params, bhcOutputs<true, true> &outputs);
#endif

}} // namespace bhc::mode
