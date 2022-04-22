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

void RayModeWorker(const bhcParams &params, bhcOutputs &outputs)
{
    rayPt<false> *localmem = nullptr;
    if(IsRayCopyMode(outputs.rayinfo)) localmem = new rayPt<false>[MaxN];
    
    try{
        
    while(true){
        int32_t job = jobID++;
        int32_t isrc, ialpha, Nsteps = -1;
        if(!GetJobIndices(isrc, ialpha, job, params.Pos, params.Angles)) break;
        if(!RunRay(outputs.rayinfo, params, localmem, job, isrc, ialpha, Nsteps)) break;
    }
    
    }catch(const std::exception &e){
        std::lock_guard<std::mutex> lock(exceptionMutex);
        exceptionStr += std::string(e.what()) + "\n";
    }
    
    if(IsRayCopyMode(outputs.rayinfo)) delete[] localmem;
}

// TL mode

cpxf *uAllSources;

void FieldModesWorker(const bhcParams &params, bhcOutputs &outputs)
{
    try{
    
    while(true){
        int32_t job = jobID++;
        int32_t isrc, ialpha;
        if(!GetJobIndices(isrc, ialpha, job, params.Pos, params.Angles)) break;
        
        real SrcDeclAngle;
        MainFieldModes(isrc, ialpha, SrcDeclAngle, outputs.uAllSources,
            params.Bdry, params.bdinfo, params.refl, params.ssp, params.Pos,
            params.Angles, params.freqinfo, params.Beam, params.beaminfo,
            outputs.eigen, outputs.arrinfo);
    }
    
    }catch(const std::exception &e){
        std::lock_guard<std::mutex> lock(exceptionMutex);
        exceptionStr += std::string(e.what()) + "\n";
    }
}

bool run_cxx(const bhcParams &params, bhcOutputs &outputs, bool singlethread)
{
    if(!api_okay) return false;
    exceptionStr = "";
    
    try{
    
    InitSelectedMode(params, outputs, singlethread);
    std::vector<std::thread> threads;
    uint32_t cores = singlethread ? 1u : bhc::max(std::thread::hardware_concurrency(), 1u);
    jobID = 0;
    for(uint32_t i=0; i<cores; ++i) threads.push_back(std::thread(
        params.Beam->RunType[0] == 'R' ? RayModeWorker : FieldModesWorker,
        std::cref(params), std::ref(outputs)));
    for(uint32_t i=0; i<cores; ++i) threads[i].join();
    
    if(!exceptionStr.empty()) throw std::runtime_error(exceptionStr);
    
    }catch(const std::exception &e){
        api_okay = false;
        PrintFileEmu &PRTFile = *(PrintFileEmu*)params.internal;
        PRTFile << "Exception caught:\n" << e.what() << "\n";
    }

    return api_okay;
}

#ifndef BHC_BUILD_CUDA
BHC_API bool run(const bhcParams &params, bhcOutputs &outputs, bool singlethread)
{
    return run_cxx(params, outputs, singlethread);
}
#endif

}
