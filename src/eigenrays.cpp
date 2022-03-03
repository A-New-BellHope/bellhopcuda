#include "eigenrays.hpp"
#include "raymode.hpp"

#include <atomic>
#include <thread>
#include <vector>

static std::atomic<int32_t> jobID;

void EigenModePostWorker(const bhcParams &params, bhcOutputs &outputs)
{
    ray2DPt *ray2D = new ray2DPt[MaxN];
    while(true){
        int32_t job = jobID++;
        if(job >= eigen->neigen) break;
        if(job >= eigen->memsize){
            printf("Had %d eigenrays but only %d fit in memory\n", eigen->neigen, eigen->memsize);
            break;
        }
        EigenHit *hit = &eigen->hits[job];
        
        memset(ray2D, 0xFE, MaxN*sizeof(ray2DPt)); //Set to garbage values for debugging
        
        real SrcDeclAngle;
        int32_t Nsteps = hit->is;
        MainRayMode(hit->isrc, hit->ialpha, SrcDeclAngle, ray2D, Nsteps,
            Bdry, bdinfo, refl, ssp, Pos, Angles, freqinfo, Beam, beaminfo);
        if(Nsteps != hit->is + 2 && Nsteps != hit->is + 3){
            printf("Eigenray isrc %d ialpha %d hit rcvr on step %d but on retrace had %d steps\n",
                hit->isrc, hit->ialpha, hit->is, Nsteps);
        }
        
        // Write the ray trajectory to RAYFile
        rayfilemutex.lock();
        WriteRay2D(SrcDeclAngle, Nsteps, RAYFile, Bdry, ray2D);
        rayfilemutex.unlock();
    }
    delete[] ray2D;
}

void FinalizeEigenMode(const bhcParams &params, bhcOutputs &outputs, 
    std::string FileRoot, bool singlethread)
{
    InitRayMode(outputs.rayinfo, params);
    
    std::cout << (int)outputs.eigen->neigen << " eigenrays\n";
    std::vector<std::thread> threads;
    uint32_t cores = singlethread ? 1u : math::max(std::thread::hardware_concurrency(), 1u);
    jobID = 0;
    for(uint32_t i=0; i<cores; ++i) threads.push_back(std::thread(EigenModePostWorker,
        params, outputs));
    for(uint32_t i=0; i<cores; ++i) threads[i].join();
    
    FinalizeRayMode(outputs.rayinfo, FileRoot, params);
}
