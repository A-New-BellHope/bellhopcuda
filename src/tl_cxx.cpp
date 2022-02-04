#include "setup.hpp"
#include "main.hpp"
#include "output.hpp"

#include <atomic>
#include <mutex>
#include <thread>
#include <vector>

std::ofstream PRTFile;
LDOFile RAYFile;
DirectOFile SHDFile;
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
EigenInfo *eigen;

static std::atomic<int32_t> jobID;


// Ray mode

static std::mutex rayfilemutex;

void RayModeWorker()
{
    ray2DPt *ray2D = new ray2DPt[MaxN];
    while(true){
        int32_t job = jobID++;
        int32_t isrc, ialpha;
        if(!GetJobIndices(isrc, ialpha, job, Pos, Angles)) break;
        
        memset(ray2D, 0xFE, MaxN*sizeof(ray2DPt)); //Set to garbage values for debugging
        
        real SrcDeclAngle;
        int32_t Nsteps = -1;
        MainRayMode(isrc, ialpha, SrcDeclAngle, ray2D, Nsteps,
            Bdry, bdinfo, refl, ssp, Pos, Angles, freqinfo, Beam, beaminfo);
        
        // Write the ray trajectory to RAYFile
        rayfilemutex.lock();
        WriteRay2D(SrcDeclAngle, Nsteps, RAYFile, Bdry, ray2D);
        rayfilemutex.unlock();
    }
    delete[] ray2D;
}

// TL mode

cpxf *uAllSources;

void FieldModesWorker()
{
    while(true){
        int32_t job = jobID++;
        int32_t isrc, ialpha;
        if(!GetJobIndices(isrc, ialpha, job, Pos, Angles)) break;
        
        real SrcDeclAngle;
        MainFieldModes(isrc, ialpha, SrcDeclAngle, uAllSources,
            Bdry, bdinfo, refl, ssp, Pos, Angles, freqinfo, Beam, beaminfo, eigen);
    }
}


// Main

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
    
    setup(FileRoot, PRTFile, RAYFile, SHDFile, Title, fT,
        Bdry, bdinfo, refl, ssp, atten, Pos, Angles, freqinfo, Beam, beaminfo, eigen);
    InitCommon(Pos, Beam);
     
    std::vector<std::thread> threads;
    uint32_t cores = singlethread ? 1u : math::max(std::thread::hardware_concurrency(), 1u);
    jobID = 0;
    
    if(Beam->RunType[0] == 'R'){
        // Ray mode
        
        for(uint32_t i=0; i<cores; ++i) threads.push_back(std::thread(RayModeWorker));
        for(uint32_t i=0; i<cores; ++i) threads[i].join();
        
    }else if(Beam->RunType[0] == 'C' || Beam->RunType[0] == 'S' || Beam->RunType[0] == 'I'){
        // TL mode
        InitTLMode(uAllSources, Pos, Beam);
        
        Stopwatch sw;
        sw.tick();
        for(uint32_t i=0; i<cores; ++i) threads.push_back(std::thread(FieldModesWorker));
        for(uint32_t i=0; i<cores; ++i) threads[i].join();
        sw.tock();
        
        FinalizeTLMode(uAllSources, SHDFile, ssp, Pos, Angles, freqinfo, Beam);
    }else if(Beam->RunType[0] == 'E'){
        // Eigenrays mode
        InitEigenMode(eigen);
        uAllSources = nullptr;
        
        Stopwatch sw;
        sw.tick();
        for(uint32_t i=0; i<cores; ++i) threads.push_back(std::thread(FieldModesWorker));
        for(uint32_t i=0; i<cores; ++i) threads[i].join();
        sw.tock();
        
        FinalizeEigenMode(Bdry, bdinfo, refl, ssp, Pos, Angles, freqinfo, Beam,
            beaminfo, eigen, RAYFile, singlethread);
    }else{
        std::cout << "Not yet implemented RunType " << Beam->RunType[0] << "\n";
        std::abort();
    }
}
