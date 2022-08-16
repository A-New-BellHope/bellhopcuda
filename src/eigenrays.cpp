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
    const bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs)
{
    rayPt<false> *localmem = nullptr;
    if(IsRayCopyMode(outputs.rayinfo)) localmem = new rayPt<false>[MaxN];
    
    try{
    
    while(true){
        uint32_t job = jobID++;
        if(job >= outputs.eigen->neigen) break;
        if(job >= outputs.eigen->memsize){
            printf("Had %d eigenrays but only %d fit in memory\n",
                outputs.eigen->neigen, outputs.eigen->memsize);
            break;
        }
        EigenHit<R3D> *hit = &outputs.eigen->hits[job];
        int32_t Nsteps = hit->is;
        RayInitInfo rinit;
        rinit.isz = hit->isz;
        rinit.ialpha = hit->ialpha;
        if constexpr(R3D){
            rinit.isx = hit->isx;
            rinit.isy = hit->isy;
            rinit.ibeta = hit->ibeta;
        }
        if(!RunRay<O3D, R3D>(outputs.rayinfo, params, localmem, job, rinit, Nsteps)){
            printf("EigenModePostWorker RunRay failed\n");
        }
        if(Nsteps != hit->is + 2 && Nsteps != hit->is + 3){
            printf("Eigenray isrc %d ialpha %d hit rcvr on step %d but on retrace had %d steps\n",
                hit->isrc, hit->ialpha, hit->is, Nsteps);
        }
    }
    
    }catch(const std::exception &e){
        std::lock_guard<std::mutex> lock(exceptionMutex);
        exceptionStr += std::string(e.what()) + "\n";
    }
    
    if(IsRayCopyMode(outputs.rayinfo)) delete[] localmem;
}

template void EigenModePostWorker<false, false>(
    const bhcParams<false, false> &params, bhcOutputs<false, false> &outputs);
template void EigenModePostWorker<true, false>(
    const bhcParams<true, false> &params, bhcOutputs<true, false> &outputs);
template void EigenModePostWorker<true, true>(
    const bhcParams<true, true> &params, bhcOutputs<true, true> &outputs);

template<bool O3D, bool R3D> void FinalizeEigenMode(
    const bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs, 
    std::string FileRoot, bool singlethread)
{
    InitRayMode<O3D, R3D>(outputs.rayinfo, params);
    
    std::cout << (int)outputs.eigen->neigen << " eigenrays\n";
    std::vector<std::thread> threads;
    uint32_t cores = singlethread ? 1u : bhc::max(std::thread::hardware_concurrency(), 1u);
    jobID = 0;
    for(uint32_t i=0; i<cores; ++i) threads.push_back(std::thread(
        EigenModePostWorker<O3D, R3D>,
        std::cref(params), std::ref(outputs)));
    for(uint32_t i=0; i<cores; ++i) threads[i].join();
    
    if(!exceptionStr.empty()) throw std::runtime_error(exceptionStr);
    
    FinalizeRayMode<O3D, R3D>(outputs.rayinfo, FileRoot, params);
}

template void FinalizeEigenMode<false, false>(
    const bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs, 
    std::string FileRoot, bool singlethread);
template void FinalizeEigenMode<true, false>(
    const bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs, 
    std::string FileRoot, bool singlethread);
template void FinalizeEigenMode<true, true>(
    const bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs, 
    std::string FileRoot, bool singlethread);

}
