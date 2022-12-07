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
#include "run.hpp"

#include <thread>
#include <vector>

namespace bhc {

std::atomic<int32_t> sharedJobID;
std::mutex exceptionMutex;
std::string exceptionStr;

template<bool O3D, bool R3D> void RayModeWorker(
    const bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs)
{
    rayPt<R3D> *localmem = nullptr;
    if(IsRayCopyMode<O3D, R3D>(outputs.rayinfo))
        localmem = (rayPt<R3D> *)malloc(MaxN * sizeof(rayPt<R3D>));

    try {
        while(true) {
            int32_t job    = sharedJobID++;
            int32_t Nsteps = -1;
            RayInitInfo rinit;
            if(!GetJobIndices<O3D>(rinit, job, params.Pos, params.Angles)) break;
            if(!RunRay<O3D, R3D>(outputs.rayinfo, params, localmem, job, rinit, Nsteps))
                break;
        }

    } catch(const std::exception &e) {
        std::lock_guard<std::mutex> lock(exceptionMutex);
        exceptionStr += std::string(e.what()) + "\n";
    }

    if(IsRayCopyMode<O3D, R3D>(outputs.rayinfo)) free(localmem);
}

#if BHC_ENABLE_2D
template void RayModeWorker<false, false>(
    const bhcParams<false, false> &params, bhcOutputs<false, false> &outputs);
#endif
#if BHC_ENABLE_NX2D
template void RayModeWorker<true, false>(
    const bhcParams<true, false> &params, bhcOutputs<true, false> &outputs);
#endif
#if BHC_ENABLE_3D
template void RayModeWorker<true, true>(
    const bhcParams<true, true> &params, bhcOutputs<true, true> &outputs);
#endif

template<bool O3D, bool R3D> inline void RunRayMode(
    bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs, uint32_t cores)
{
    GlobalLog("%d threads\n", cores);
    std::vector<std::thread> threads;
    for(uint32_t i = 0; i < cores; ++i)
        threads.push_back(
            std::thread(RayModeWorker<O3D, R3D>, std::ref(params), std::ref(outputs)));
    for(uint32_t i = 0; i < cores; ++i) threads[i].join();
}

template<bool O3D, bool R3D> bool run(
    bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs, bool singlethread)
{
    if(!api_okay) return false;
    exceptionStr = "";
    sharedJobID  = 0;

#ifdef BHC_BUILD_CUDA
    if(singlethread) {
        GlobalLog("Single threaded mode is nonsense on CUDA, ignoring");
        singlethread = false;
    }
#endif
    uint32_t cores = std::thread::hardware_concurrency();
    if(singlethread || cores < 1u) cores = 1u;

    try {
        InitSelectedMode<O3D, R3D>(params, outputs, singlethread);
        if(IsRayRun(params.Beam)) {
            RunRayMode(params, outputs, cores);
        } else {
            RunFieldModesSelInfl(params, outputs, cores);
        }
        if(!exceptionStr.empty()) throw std::runtime_error(exceptionStr);
    } catch(const std::exception &e) {
        api_okay              = false;
        PrintFileEmu &PRTFile = *(PrintFileEmu *)params.internal;
        PRTFile << "Exception caught:\n" << e.what() << "\n";
    }

    return api_okay;
}

#if BHC_ENABLE_2D
template bool BHC_API run<false, false>(
    bhcParams<false, false> &params, bhcOutputs<false, false> &outputs,
    bool singlethread);
#endif
#if BHC_ENABLE_NX2D
template bool BHC_API run<true, false>(
    bhcParams<true, false> &params, bhcOutputs<true, false> &outputs, bool singlethread);
#endif
#if BHC_ENABLE_3D
template bool BHC_API run<true, true>(
    bhcParams<true, true> &params, bhcOutputs<true, true> &outputs, bool singlethread);
#endif

} // namespace bhc
