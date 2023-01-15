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
    const bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs, ErrState *errState)
{
    SetupThread();

    rayPt<R3D> *localmem = nullptr;
    if(IsRayCopyMode<O3D, R3D>(outputs.rayinfo))
        localmem = (rayPt<R3D> *)malloc(MaxN * sizeof(rayPt<R3D>));

    try {
        while(true) {
            int32_t job    = GetInternal(params)->sharedJobID++;
            int32_t Nsteps = -1;
            RayInitInfo rinit;
            if(!GetJobIndices<O3D>(rinit, job, params.Pos, params.Angles)) break;
            if(!RunRay<O3D, R3D>(
                   outputs.rayinfo, params, localmem, job, rinit, Nsteps, errState)) {
                break;
            }
        }

    } catch(const std::exception &e) {
        std::lock_guard<std::mutex> lock(GetInternal(params)->exceptionMutex);
        GetInternal(params)->exceptionStr += std::string(e.what()) + "\n";
    }

    if(IsRayCopyMode<O3D, R3D>(outputs.rayinfo)) free(localmem);
}

#if BHC_ENABLE_2D
template void RayModeWorker<false, false>(
    const bhcParams<false, false> &params, bhcOutputs<false, false> &outputs,
    ErrState *errState);
#endif
#if BHC_ENABLE_NX2D
template void RayModeWorker<true, false>(
    const bhcParams<true, false> &params, bhcOutputs<true, false> &outputs,
    ErrState *errState);
#endif
#if BHC_ENABLE_3D
template void RayModeWorker<true, true>(
    const bhcParams<true, true> &params, bhcOutputs<true, true> &outputs,
    ErrState *errState);
#endif

template<bool O3D, bool R3D> inline void RunRayMode(
    bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs)
{
    ErrState errState;
    ResetErrState(&errState);
    uint32_t nthreads = GetNumThreads(params.maxThreads);
    EXTWARN(
        "%d threads, copy mode %s", nthreads,
        IsRayCopyMode<O3D, R3D>(outputs.rayinfo) ? "true" : "false");
    std::vector<std::thread> threads;
    for(uint32_t i = 0; i < nthreads; ++i)
        threads.push_back(std::thread(
            RayModeWorker<O3D, R3D>, std::ref(params), std::ref(outputs), &errState));
    for(uint32_t i = 0; i < nthreads; ++i) threads[i].join();
    CheckReportErrors(GetInternal(params), &errState);
}

template<bool O3D, bool R3D> bool run(
    bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs)
{
    if(!GetInternal(params)->api_okay) return false;
    GetInternal(params)->exceptionStr = "";
    GetInternal(params)->sharedJobID  = 0;

    try {
        Stopwatch sw;
        sw.tick();
        InitSelectedMode<O3D, R3D>(params, outputs);
        sw.tock("InitSelectedMode");
        if(IsRayRun(params.Beam)) {
            RunRayMode(params, outputs);
        } else {
#ifdef BHC_BUILD_CUDA
            checkCudaErrors(cudaSetDevice(GetInternal(params)->m_gpu));
#endif
            RunFieldModesSelInfl(params, outputs);
        }
        if(!GetInternal(params)->exceptionStr.empty()) {
            throw std::runtime_error(GetInternal(params)->exceptionStr);
        }
    } catch(const std::exception &e) {
        GetInternal(params)->api_okay = false;
        PrintFileEmu &PRTFile         = GetInternal(params)->PRTFile;
        PRTFile << "Exception caught:\n" << e.what() << "\n";
    }

    return GetInternal(params)->api_okay;
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
    const bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs)
{
    if(!GetInternal(params)->api_okay) return false;

    if(IsRayRun(params.Beam)) {
        // Ray mode
        bhc::WriteOutRays<O3D, R3D>(params, outputs.rayinfo);
    } else if(IsTLRun(params.Beam)) {
        // TL mode
        bhc::WriteOutTL<O3D, R3D>(params, outputs);
    } else if(IsEigenraysRun(params.Beam)) {
        // Eigenrays mode
        bhc::WriteOutEigenrays<O3D, R3D>(params, outputs);
    } else if(IsArrivalsRun(params.Beam)) {
        // Arrivals mode
        bhc::WriteOutArrivals<O3D, R3D>(
            params.Pos, params.freqinfo, params.Beam, GetInternal(params)->FileRoot,
            outputs.arrinfo);
    } else {
        std::cout << "Invalid RunType " << params.Beam->RunType[0] << "\n";
        GetInternal(params)->api_okay = false;
    }

    return GetInternal(params)->api_okay;
}

#if BHC_ENABLE_2D
template BHC_API bool writeout<false, false>(
    const bhcParams<false, false> &params, bhcOutputs<false, false> &outputs);
#if BHC_ENABLE_NX2D
#endif
template BHC_API bool writeout<true, false>(
    const bhcParams<true, false> &params, bhcOutputs<true, false> &outputs);
#if BHC_ENABLE_3D
#endif
template BHC_API bool writeout<true, true>(
    const bhcParams<true, true> &params, bhcOutputs<true, true> &outputs);
#endif

} // namespace bhc
