#include "eigenrays.hpp"
#include "step.hpp"
#include "main.hpp"
#include "output.hpp"

#include <atomic>
#include <mutex>
#include <thread>
#include <vector>

static std::atomic<int32_t> jobID;

static std::mutex rayfilemutex;

void EigenModePostWorker(
    const BdryType *Bdry, const BdryInfo *bdinfo, const ReflectionInfo *refl,
    const SSPStructure *ssp, const Position *Pos, const AnglesStructure *Angles,
    const FreqInfo *freqinfo, const BeamStructure *Beam, const BeamInfo *beaminfo,
    const EigenInfo *eigen, LDOFile &RAYFile)
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

void FinalizeEigenMode(
    const BdryType *Bdry, const BdryInfo *bdinfo, const ReflectionInfo *refl,
    const SSPStructure *ssp, const Position *Pos, const AnglesStructure *Angles,
    const FreqInfo *freqinfo, const BeamStructure *Beam, const BeamInfo *beaminfo,
    const EigenInfo *eigen, LDOFile &RAYFile, bool singlethread)
{
    std::cout << (int)eigen->neigen << " eigenrays\n";
    std::vector<std::thread> threads;
    uint32_t cores = singlethread ? 1u : math::max(std::thread::hardware_concurrency(), 1u);
    jobID = 0;
    for(uint32_t i=0; i<cores; ++i) threads.push_back(std::thread(EigenModePostWorker,
        Bdry, bdinfo, refl, ssp, Pos, Angles, freqinfo, Beam, beaminfo, eigen,
        std::ref(RAYFile)));
    for(uint32_t i=0; i<cores; ++i) threads[i].join();
}
