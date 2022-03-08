#include "run.hpp"

#include <atomic>
#include <thread>
#include <vector>

namespace bhc {

static std::atomic<int32_t> jobID;

// Ray mode

void RayModeWorker(const bhcParams &params, bhcOutputs &outputs)
{
    ray2DPt *localmem = nullptr;
    if(IsRayCopyMode(outputs.rayinfo)) localmem = new ray2DPt[MaxN];
    while(true){
        int32_t job = jobID++;
        int32_t isrc, ialpha, Nsteps = -1;
        if(!GetJobIndices(isrc, ialpha, job, params.Pos, params.Angles)) break;
        if(!RunRay(outputs.rayinfo, params, localmem, job, isrc, ialpha, Nsteps)) break;
    }
    if(IsRayCopyMode(outputs.rayinfo)) delete[] localmem;
}

// TL mode

cpxf *uAllSources;

void FieldModesWorker(const bhcParams &params, bhcOutputs &outputs)
{
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
}

void run_cxx(std::ostream &PRTFile, const bhcParams &params, bhcOutputs &outputs,
    bool singlethread)
{
    InitSelectedMode(PRTFile, params, outputs, singlethread);
    std::vector<std::thread> threads;
    uint32_t cores = singlethread ? 1u : bhc::max(std::thread::hardware_concurrency(), 1u);
    jobID = 0;
    for(uint32_t i=0; i<cores; ++i) threads.push_back(std::thread(
        params.Beam->RunType[0] == 'R' ? RayModeWorker : FieldModesWorker,
        std::cref(params), std::ref(outputs)));
    for(uint32_t i=0; i<cores; ++i) threads[i].join();
}

#ifndef BHC_BUILD_CUDA
BHC_API void run(std::ostream &PRTFile, const bhcParams &params, bhcOutputs &outputs,
    bool singlethread)
{
    run_cxx(PRTFile, params, outputs, singlethread);
}
#endif

}
