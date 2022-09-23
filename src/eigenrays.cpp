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

#include <atomic>
#include <mutex>
#include <thread>
#include <vector>

namespace bhc {

static std::atomic<uint32_t> jobID;
static std::mutex exceptionMutex;
static std::string exceptionStr;

template<bool O3D, bool R3D> void EigenModePostWorker(
    bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs)
{
    rayPt<R3D> *localmem = nullptr;
    if(IsRayCopyMode(outputs.rayinfo)) localmem = new rayPt<R3D>[MaxN];
    
    try{
    
    while(true){
        uint32_t job = jobID++;
        if(job >= outputs.eigen->neigen) break;
        if(job >= outputs.eigen->memsize){
            GlobalLog("Had %d eigenrays but only %d fit in memory\n",
                outputs.eigen->neigen, outputs.eigen->memsize);
            break;
        }
        EigenHit *hit = &outputs.eigen->hits[job];
        int32_t Nsteps = hit->is;
        RayInitInfo rinit;
        rinit.isx = hit->isx;
        rinit.isy = hit->isy;
        rinit.isz = hit->isz;
        rinit.ialpha = hit->ialpha;
        rinit.ibeta = hit->ibeta;
        if(!RunRay<O3D, R3D>(outputs.rayinfo, params, localmem, job, rinit, Nsteps)){
            GlobalLog("EigenModePostWorker RunRay failed\n");
        }
        if(Nsteps != hit->is + 2 && Nsteps != hit->is + 3){
            GlobalLog("Eigenray isxyz (%d,%d,%d) ialpha/beta (%d,%d) "
                "hit rcvr on step %d but on retrace had %d steps\n",
                hit->isx, hit->isy, hit->isz, hit->ialpha, hit->ibeta, hit->is, Nsteps);
        }
    }
    
    }catch(const std::exception &e){
        std::lock_guard<std::mutex> lock(exceptionMutex);
        exceptionStr += std::string(e.what()) + "\n";
    }
    
    if(IsRayCopyMode(outputs.rayinfo)) delete[] localmem;
}

template void EigenModePostWorker<false, false>(
    bhcParams<false, false> &params, bhcOutputs<false, false> &outputs);
template void EigenModePostWorker<true, false>(
    bhcParams<true, false> &params, bhcOutputs<true, false> &outputs);
template void EigenModePostWorker<true, true>(
    bhcParams<true, true> &params, bhcOutputs<true, true> &outputs);

template<bool O3D, bool R3D> void FinalizeEigenMode(
    bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs, 
    std::string FileRoot, bool singlethread)
{
    InitRayMode<O3D, R3D>(outputs.rayinfo, params);
    
    GlobalLog("%d eigenrays\n", (int)outputs.eigen->neigen);
    std::vector<std::thread> threads;
    uint32_t cores = singlethread ? 1u : bhc::max(std::thread::hardware_concurrency(), 1u);
    jobID = 0;
    for(uint32_t i=0; i<cores; ++i) threads.push_back(std::thread(
        EigenModePostWorker<O3D, R3D>,
        std::ref(params), std::ref(outputs)));
    for(uint32_t i=0; i<cores; ++i) threads[i].join();
    
    if(!exceptionStr.empty()) throw std::runtime_error(exceptionStr);
    
    FinalizeRayMode<O3D, R3D>(outputs.rayinfo, FileRoot, params);
}

template void FinalizeEigenMode<false, false>(
    bhcParams<false, false> &params, bhcOutputs<false, false> &outputs, 
    std::string FileRoot, bool singlethread);
template void FinalizeEigenMode<true, false>(
    bhcParams<true, false> &params, bhcOutputs<true, false> &outputs, 
    std::string FileRoot, bool singlethread);
template void FinalizeEigenMode<true, true>(
    bhcParams<true, true> &params, bhcOutputs<true, true> &outputs, 
    std::string FileRoot, bool singlethread);

}
