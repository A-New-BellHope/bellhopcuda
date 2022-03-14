#include "eigenrays.hpp"
#include "raymode.hpp"

#include <atomic>
#include <thread>
#include <vector>

namespace bhc {

static std::atomic<uint32_t> jobID;

void EigenModePostWorker(const bhcParams &params, bhcOutputs &outputs)
{
    ray2DPt *localmem = nullptr;
    if(IsRayCopyMode(outputs.rayinfo)) localmem = new ray2DPt[MaxN];
    while(true){
        uint32_t job = jobID++;
        if(job >= outputs.eigen->neigen) break;
        if(job >= outputs.eigen->memsize){
            printf("Had %d eigenrays but only %d fit in memory\n",
                outputs.eigen->neigen, outputs.eigen->memsize);
            break;
        }
        EigenHit *hit = &outputs.eigen->hits[job];
        int32_t Nsteps = hit->is;
        if(!RunRay(outputs.rayinfo, params, localmem, job, hit->isrc, hit->ialpha, Nsteps)){
            printf("EigenModePostWorker RunRay failed\n");
        }
        if(Nsteps != hit->is + 2 && Nsteps != hit->is + 3){
            printf("Eigenray isrc %d ialpha %d hit rcvr on step %d but on retrace had %d steps\n",
                hit->isrc, hit->ialpha, hit->is, Nsteps);
        }
    }
    if(IsRayCopyMode(outputs.rayinfo)) delete[] localmem;
}

void FinalizeEigenMode(const bhcParams &params, bhcOutputs &outputs, 
    std::string FileRoot, bool singlethread)
{
    InitRayMode(outputs.rayinfo, params);
    
    std::cout << (int)outputs.eigen->neigen << " eigenrays\n";
    std::vector<std::thread> threads;
    uint32_t cores = singlethread ? 1u : bhc::max(std::thread::hardware_concurrency(), 1u);
    jobID = 0;
    for(uint32_t i=0; i<cores; ++i) threads.push_back(std::thread(EigenModePostWorker,
        std::cref(params), std::ref(outputs)));
    for(uint32_t i=0; i<cores; ++i) threads[i].join();
    
    FinalizeRayMode(outputs.rayinfo, FileRoot, params);
}

}
