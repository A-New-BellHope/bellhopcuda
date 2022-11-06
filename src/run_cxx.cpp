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

#include <atomic>
#include <mutex>
#include <thread>
#include <vector>

namespace bhc {

static std::atomic<int32_t> jobID;
static std::mutex exceptionMutex;
static std::string exceptionStr;

// Ray mode

template<bool O3D, bool R3D> void RayModeWorker(
    const bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs)
{
    rayPt<R3D> *localmem = nullptr;
    if(IsRayCopyMode<O3D, R3D>(outputs.rayinfo)) localmem = new rayPt<R3D>[MaxN];

    try {
        while(true) {
            int32_t job    = jobID++;
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

    if(IsRayCopyMode<O3D, R3D>(outputs.rayinfo)) delete[] localmem;
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

// TL, arrivals, eigenrays

cpxf *uAllSources;

// LP: Want this one to be where the different functions are instantiated, and
// the functions below to be inlined.
template<typename RI, bool O3D, bool R3D> CPU_NOINLINE void FieldModesWorkerInner(
    const bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs)
{
    while(true) {
        int32_t job = jobID++;
        RayInitInfo rinit;
        if(!GetJobIndices<O3D>(rinit, job, params.Pos, params.Angles)) break;

        MainFieldModes<RI, O3D, R3D>(
            rinit, outputs.uAllSources, params.Bdry, params.bdinfo, params.refl,
            params.ssp, params.Pos, params.Angles, params.freqinfo, params.Beam,
            params.beaminfo, outputs.eigen, outputs.arrinfo);
    }
}

template<char IT, bool O3D, bool R3D> inline void RunSelectRunType(
    const bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs)
{
    char rt = params.Beam->RunType[0];
    if(rt == 'C') {
        FieldModesWorkerInner<RIType<'C', IT>, O3D, R3D>(params, outputs);
    } else if(rt == 'S') {
        FieldModesWorkerInner<RIType<'S', IT>, O3D, R3D>(params, outputs);
    } else if(rt == 'I') {
        FieldModesWorkerInner<RIType<'I', IT>, O3D, R3D>(params, outputs);
    } else if(rt == 'E') {
        FieldModesWorkerInner<RIType<'E', IT>, O3D, R3D>(params, outputs);
    } else if(rt == 'A') {
        FieldModesWorkerInner<RIType<'A', IT>, O3D, R3D>(params, outputs);
    } else if(rt == 'a') {
        FieldModesWorkerInner<RIType<'a', IT>, O3D, R3D>(params, outputs);
    } else {
        GlobalLog("Invalid RunType for FieldModesWorker %c!", rt);
        bail();
    }
}

template<bool O3D, bool R3D> inline void RunSelectInflType(
    const bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs)
{
    char it = params.Beam->Type[0];
    if(it == 'R') {
        RunSelectRunType<'R', O3D, R3D>(params, outputs);
    } else if(it == 'C') {
        RunSelectRunType<'C', O3D, R3D>(params, outputs);
    } else if(it == 'G' || it == '^' || it == ' ') {
        RunSelectRunType<'G', O3D, R3D>(params, outputs);
    } else if(it == 'g') {
        RunSelectRunType<'g', O3D, R3D>(params, outputs);
    } else if(it == 'B') {
        RunSelectRunType<'B', O3D, R3D>(params, outputs);
    } else if(it == 'b') {
        RunSelectRunType<'b', O3D, R3D>(params, outputs);
    } else if(it == 'S') {
        RunSelectRunType<'S', O3D, R3D>(params, outputs);
    } else {
        GlobalLog("Invalid Beam type %c!", it);
        bail();
    }
}

template<bool O3D, bool R3D> CPU_NOINLINE void FieldModesWorker(
    const bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs)
{
    try {
        RunSelectInflType<O3D, R3D>(params, outputs);
    } catch(const std::exception &e) {
        std::lock_guard<std::mutex> lock(exceptionMutex);
        exceptionStr += std::string(e.what()) + "\n";
    }
}

#if BHC_ENABLE_2D
template void FieldModesWorker<false, false>(
    const bhcParams<false, false> &params, bhcOutputs<false, false> &outputs);
#endif
#if BHC_ENABLE_NX2D
template void FieldModesWorker<true, false>(
    const bhcParams<true, false> &params, bhcOutputs<true, false> &outputs);
#endif
#if BHC_ENABLE_3D
template void FieldModesWorker<true, true>(
    const bhcParams<true, true> &params, bhcOutputs<true, true> &outputs);
#endif

template<bool O3D, bool R3D> bool run_cxx(
    bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs, bool singlethread)
{
    if(!api_okay) return false;
    exceptionStr = "";

    try {
        InitSelectedMode<O3D, R3D>(params, outputs, singlethread);
        std::vector<std::thread> threads;
        uint32_t cores = singlethread ? 1u
                                      : bhc::max(std::thread::hardware_concurrency(), 1u);
        jobID          = 0;
        for(uint32_t i = 0; i < cores; ++i)
            threads.push_back(std::thread(
                IsRayRun(params.Beam) ? RayModeWorker<O3D, R3D>
                                      : FieldModesWorker<O3D, R3D>,
                std::ref(params), std::ref(outputs)));
        for(uint32_t i = 0; i < cores; ++i) threads[i].join();

        if(!exceptionStr.empty()) throw std::runtime_error(exceptionStr);

    } catch(const std::exception &e) {
        api_okay              = false;
        PrintFileEmu &PRTFile = *(PrintFileEmu *)params.internal;
        PRTFile << "Exception caught:\n" << e.what() << "\n";
    }

    return api_okay;
}

#if BHC_ENABLE_2D
template bool run_cxx<false, false>(
    bhcParams<false, false> &params, bhcOutputs<false, false> &outputs,
    bool singlethread);
#endif
#if BHC_ENABLE_NX2D
template bool run_cxx<true, false>(
    bhcParams<true, false> &params, bhcOutputs<true, false> &outputs, bool singlethread);
#endif
#if BHC_ENABLE_3D
template bool run_cxx<true, true>(
    bhcParams<true, true> &params, bhcOutputs<true, true> &outputs, bool singlethread);
#endif

#ifndef BHC_BUILD_CUDA

template<bool O3D, bool R3D> bool run(
    bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs, bool singlethread)
{
    return run_cxx(params, outputs, singlethread);
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

#endif

} // namespace bhc
