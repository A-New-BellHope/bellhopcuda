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
#include "eigenrays.hpp"
#include "raymode.hpp"
#include "run.hpp"

#include <vector>

namespace bhc {

template<bool O3D, bool R3D> void EigenModePostWorker(
    bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs)
{
    SetupThread();

    rayPt<R3D> *localmem = nullptr;
    if(IsRayCopyMode<O3D, R3D>(outputs.rayinfo))
        localmem = (rayPt<R3D> *)malloc(MaxN * sizeof(rayPt<R3D>));

    try {
        while(true) {
            uint32_t job = sharedJobID++;
            if(job >= outputs.eigen->neigen) break;
            if(job >= outputs.eigen->memsize) {
                GlobalLog(
                    "Had %d eigenrays but only %d fit in memory\n", outputs.eigen->neigen,
                    outputs.eigen->memsize);
                break;
            }
            EigenHit *hit  = &outputs.eigen->hits[job];
            int32_t Nsteps = hit->is;
            RayInitInfo rinit;
            rinit.isx    = hit->isx;
            rinit.isy    = hit->isy;
            rinit.isz    = hit->isz;
            rinit.ialpha = hit->ialpha;
            rinit.ibeta  = hit->ibeta;
            if(!RunRay<O3D, R3D>(outputs.rayinfo, params, localmem, job, rinit, Nsteps)) {
                GlobalLog("EigenModePostWorker RunRay failed\n");
            }
        }

    } catch(const std::exception &e) {
        std::lock_guard<std::mutex> lock(exceptionMutex);
        exceptionStr += std::string(e.what()) + "\n";
    }

    if(IsRayCopyMode<O3D, R3D>(outputs.rayinfo)) free(localmem);
}

#if BHC_ENABLE_2D
template void EigenModePostWorker<false, false>(
    bhcParams<false, false> &params, bhcOutputs<false, false> &outputs);
#endif
#if BHC_ENABLE_NX2D
template void EigenModePostWorker<true, false>(
    bhcParams<true, false> &params, bhcOutputs<true, false> &outputs);
#endif
#if BHC_ENABLE_3D
template void EigenModePostWorker<true, true>(
    bhcParams<true, true> &params, bhcOutputs<true, true> &outputs);
#endif

template<bool O3D, bool R3D> void WriteOutEigenrays(
    bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs, std::string FileRoot)
{
    InitRayMode<O3D, R3D>(outputs.rayinfo, params, outputs.eigen->neigen);

    GlobalLog("%d eigenrays\n", (int)outputs.eigen->neigen);
    std::vector<std::thread> threads;
    sharedJobID       = 0;
    uint32_t nthreads = GetNumThreads(params.maxThreads);
    for(uint32_t i = 0; i < nthreads; ++i)
        threads.push_back(std::thread(
            EigenModePostWorker<O3D, R3D>, std::ref(params), std::ref(outputs)));
    for(uint32_t i = 0; i < nthreads; ++i) threads[i].join();

    if(!exceptionStr.empty()) throw std::runtime_error(exceptionStr);

    WriteOutRays<O3D, R3D>(outputs.rayinfo, FileRoot, params);
}

#if BHC_ENABLE_2D
template void WriteOutEigenrays<false, false>(
    bhcParams<false, false> &params, bhcOutputs<false, false> &outputs,
    std::string FileRoot);
#endif
#if BHC_ENABLE_NX2D
template void WriteOutEigenrays<true, false>(
    bhcParams<true, false> &params, bhcOutputs<true, false> &outputs,
    std::string FileRoot);
#endif
#if BHC_ENABLE_3D
template void WriteOutEigenrays<true, true>(
    bhcParams<true, true> &params, bhcOutputs<true, true> &outputs, std::string FileRoot);
#endif

} // namespace bhc
