#include "setup.hpp"
#include "trace.hpp"
#include "output.hpp"

#include <atomic>
#include <mutex>
#include <thread>
#include <vector>

std::ofstream PRTFile, ARRFile;
LDOFile RAYFile;
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
        int32_t isrc, ialpha; // LP: `is` changed to `isrc` because `is` is used for steps
        if(Angles->iSingle_alpha >= 0){
            isrc = ray;
            ialpha = Angles->iSingle_alpha;
        }else{
            isrc = ray / Angles->Nalpha;
            ialpha = ray % Angles->Nalpha;
        }
        if(isrc >= Pos->NSz) break;
        
        memset(ray2D, 0xFE, MaxN*sizeof(ray2DPt)); //Set to garbage values for debugging
        
        real SrcDeclAngle;
        int32_t Nsteps;
        MainRayMode(isrc, ialpha, SrcDeclAngle, ray2D, Nsteps, 
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
    std::string FileRoot;
    bool singlethread = false;
    for(int32_t i=1; i<argc; ++i){
        std::string s = argv[i];
        if(argv[i][0] == '-'){
            if(s.length() >= 2 && argv[i][1] == '-'){ //two dashes
                s = s.substr(1);
            }
            if(s == "-1" || s == "-singlethread"){
                singlethread = true;
            }else{
                std::cout << "Unknown command-line option \"" << s << "\"\n";
                std::abort();
            }
        }else{
            if(FileRoot.empty()){
                FileRoot = s;
            }else{
                std::cout << "Intepreting both \"" << FileRoot << "\" and \"" << 
                    s << "\" as FileRoot, error\n";
                std::abort();
            }
        }
    }
    if(FileRoot.empty()){
        std::cout << "Must provide FileRoot as command-line parameter\n";
        std::abort();
    }
    
    setup(FileRoot, PRTFile, RAYFile, ARRFile, Title, fT,
        Bdry, bdinfo, refl, ssp, atten, Pos, Angles, freqinfo, Beam, beaminfo);
    core_setup(PRTFile, fT, Bdry, bdinfo, atten, Angles, freqinfo, Beam);
    
    rayID = 0;
    
    std::vector<std::thread> threads;
    uint32_t cores = singlethread ? 1u : std::max(std::thread::hardware_concurrency(), 1u);
    for(uint32_t i=0; i<cores; ++i) threads.push_back(std::thread(RayWorker));
    for(uint32_t i=0; i<cores; ++i) threads[i].join();
}
