#include "setup.hpp"
#include "trace.hpp"

#include <atomic>
#include <mutex>
#include <vector>

std::ofstream PRTFile, RAYFile, ARRFile;
std::string Title;
real fT;
BdryType *Bdry;
BdryInfo *bdinfo;
ReflectionInfo *refl;
SSPStructure *ssp;
AttenInfo *atten;
Position *Pos;
AnglesStructure *Angles;
FreqInfo *freqinfo;
BeamStructure *Beam;
BeamInfo *beaminfo;

std::atomic<int32_t> rayID;
std::mutex rayfilemutex;

void RayWorker()
{
    ray2DPt *ray2D = new ray2DPt[MaxN];
    while(true){
        int32_t ray = rayID++;
        int32_t is, ialpha;
        if(Angles->iSingleAlpha >= 0){
            is = ray;
            ialpha = Angles->iSingleAlpha;
        }else{
            is = ray / Angles->Nalpha;
            ialpha = ray % Angles->Nalpha;
        }
        if(is >= Pos->NSz) break;
        
        memset(ray2D, 0xFE, MaxN*sizeof(ray2DPt)); //Set to garbage values for debugging
        
        real SrcDeclAngle;
        int32_t Nsteps;
        CoreSingleBeam(is, ialpha, ray2D, SrcDeclAngle, Nsteps, 
            Bdry, bdinfo, refl, ssp, Pos, Angles, freqinfo, Beam, beaminfo);
        
        if(Beam->RunType[0] == 'R'){
            // Write the ray trajectory to RAYFile
            rayfilemutex.lock();
            WriteRay2D(SrcDeclAngle, Nsteps, RAYFile, Bdry, ray2D);
            rayfilemutex.unlock();
        }else{
            // Compute the contribution to the field
            std::cout << "TODO Influence not yet implemented\n";
            std::abort();
        }
    }
    delete[] ray2D;
}

int main(int argc, char **argv)
{
    
    setup(argc, argv, PRTFile, RAYFile, ARRFile, Title, fT,
        Bdry, bdinfo, refl, ssp, atten, Pos, Angles, freqinfo, Beam, beaminfo);
    core_setup(PRTFile, Bdry, bdinfo, atten, Angles, freqinfo, Beam);
    
    rayID = 0;
    
    std::vector<std::thread> threads;
    int cores = std::max(std::thread::hardware_concurrency(), 1);
    for(int i=0; i<cores; ++i) threads.push_back(std::thread(RayWorker));
    for(int i=0; i<cores; ++i) threads[i].join();
}
