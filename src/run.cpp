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
#include "run.hpp"

#include <vector>

namespace bhc {

template<bool O3D, bool R3D> void RayModeWorker(
    const bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs, int32_t worker,
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
    const bhcParams<false, false> &params, bhcOutputs<false, false> &outputs,
    int32_t worker, ErrState *errState);
#endif
#if BHC_ENABLE_NX2D
template void RayModeWorker<true, false>(
    const bhcParams<true, false> &params, bhcOutputs<true, false> &outputs,
    int32_t worker, ErrState *errState);
#endif
#if BHC_ENABLE_3D
template void RayModeWorker<true, true>(
    const bhcParams<true, true> &params, bhcOutputs<true, true> &outputs, int32_t worker,
    ErrState *errState);
#endif

template<bool O3D, bool R3D> inline void RunRayMode(
    bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs)
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

template<bool O3D, bool R3D> bool run(
    bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs)
{
    try {
        Stopwatch sw(GetInternal(params));
        sw.tick();
        InitSelectedMode<O3D, R3D>(params, outputs);
        sw.tock("InitSelectedMode");

        sw.tick();
        if(IsRayRun(params.Beam)) {
            RunRayMode(params, outputs);
            sw.tock("RunRayMode");
        } else {
#ifdef BHC_BUILD_CUDA
            checkCudaErrors(cudaSetDevice(GetInternal(params)->gpuIndex));
#endif
            RunFieldModesSelInfl(params, outputs);
            sw.tock("RunFieldModes");
        }

        sw.tick();
        if(IsRayRun(params.Beam)) {
            bhc::PostProcessRays<O3D, R3D>(params, outputs.rayinfo);
        } else if(IsTLRun(params.Beam)) {
            bhc::PostProcessTL<O3D, R3D>(params, outputs);
        } else if(IsEigenraysRun(params.Beam)) {
            bhc::PostProcessEigenrays<O3D, R3D>(params, outputs);
        } else if(IsArrivalsRun(params.Beam)) {
            bhc::PostProcessArrivals<O3D, R3D>(params, outputs.arrinfo);
        } else {
            EXTERR("Invalid RunType %c\n", params.Beam->RunType[0]);
        }
        sw.tock("PostProcess");

    } catch(const std::exception &e) {
        EXTWARN("Exception caught in bhc::run(): %s\n", e.what());
        return false;
    }

    return true;
}

#if BHC_ENABLE_2D
template bool BHC_API
run<false, false>(bhcParams<false, false> &params, bhcOutputs<false, false> &outputs);
#endif
#if BHC_ENABLE_NX2D
template bool BHC_API
run<true, false>(bhcParams<true, false> &params, bhcOutputs<true, false> &outputs);
#endif
#if BHC_ENABLE_3D
template bool BHC_API
run<true, true>(bhcParams<true, true> &params, bhcOutputs<true, true> &outputs);
#endif

template<bool O3D, bool R3D> bool writeout(
    const bhcParams<O3D, R3D> &params, const bhcOutputs<O3D, R3D> &outputs)
{
    try {
        Stopwatch sw(GetInternal(params));
        sw.tick();
        if(IsRayRun(params.Beam)) {
            bhc::WriteOutRays<O3D, R3D>(params, outputs.rayinfo);
        } else if(IsTLRun(params.Beam)) {
            bhc::WriteOutTL<O3D, R3D>(params, outputs);
        } else if(IsEigenraysRun(params.Beam)) {
            bhc::WriteOutEigenrays<O3D, R3D>(params, outputs);
        } else if(IsArrivalsRun(params.Beam)) {
            bhc::WriteOutArrivals<O3D, R3D>(params, outputs.arrinfo);
        } else {
            EXTERR("Invalid RunType %c\n", params.Beam->RunType[0]);
        }
        sw.tock("writeout");
    } catch(const std::exception &e) {
        EXTWARN("Exception caught in bhc::writeout(): %s\n", e.what());
        return false;
    }

    return true;
}

#if BHC_ENABLE_2D
template BHC_API bool writeout<false, false>(
    const bhcParams<false, false> &params, const bhcOutputs<false, false> &outputs);
#endif
#if BHC_ENABLE_NX2D
template BHC_API bool writeout<true, false>(
    const bhcParams<true, false> &params, const bhcOutputs<true, false> &outputs);
#endif
#if BHC_ENABLE_3D
template BHC_API bool writeout<true, true>(
    const bhcParams<true, true> &params, const bhcOutputs<true, true> &outputs);
#endif

} // namespace bhc
