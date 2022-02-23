//Joe Snider
//12/21
//
//Testing dll target for the bellhop translation.
//Not real general.

#define BELLHOPCXX_EXPORTS

#include "setup.hpp"
#include "trace.hpp"
#include "output.hpp"
#include "trace_dll.hpp"

#include <atomic>
#include <mutex>
#include <thread>
#include <vector>
#include <fstream>

std::ofstream PRTFile, ARRFile;
LDOFile RAYFile;
std::string Title;
real fT;
BdryType* Bdry;
BdryInfo* bdinfo;
ReflectionInfo* refl;
SSPStructure* ssp;
AttenInfo* atten;
Position* Pos;
AnglesStructure* Angles;
FreqInfo* freqinfo;
BeamStructure* Beam;
BeamInfo* beaminfo;

std::atomic<int32_t> rayID;
std::mutex rayfilemutex;


//just dump them
std::ostream& operator<<(std::ostream& out, const ray2DPt& x) {
    out << x.NumTopBnc << " "
        << x.NumBotBnc << " "
        << x.x.x << " " << x.x.y << " "
        << x.t.x << " " << x.t.y << " "
        << x.p.x << " " << x.p.y << " "
        << x.q.x << " " << x.q.y << " "
        << x.c << " "
        << x.Amp << " " << x.Phase << " "
        << x.tau;
    return out;
}

//thread communication lazyness
//make sure there's enough room for MaxN * Pos->NSz values
ray2DPt *ray2DAll;

void RayWorker()
{
    while (true) {
        int32_t ray = rayID;
        int32_t is, ialpha;
        ray2DPt* ray2D = ray2DAll + MaxN * ray;
        rayID += 1;
        //{
        //    std::lock_guard<std::mutex> g(rayfilemutex);
        //    std::cout << "Starting thread " << ray 
        //        << " at " << reinterpret_cast<void*>(ray2D) 
        //        << "\n" << std::flush;
        //}
        if (Angles->iSingle_alpha >= 0) {
            is = ray;
            ialpha = Angles->iSingle_alpha;
        }
        else {
            is = ray / Angles->Nalpha;
            ialpha = ray % Angles->Nalpha;
        }
        if (is >= Pos->NSz) break;

        memset(ray2D, 0xFE, MaxN * sizeof(ray2DPt)); //Set to garbage values for debugging

        real SrcDeclAngle;
        int32_t Nsteps;
        CoreSingleBeam(is, ialpha, ray2D, SrcDeclAngle, Nsteps,
            Bdry, bdinfo, refl, ssp, Pos, Angles, freqinfo, Beam, beaminfo);

        //{
        //    std::lock_guard<std::mutex> g(rayfilemutex);
        //    std::cout << "Ending thread " << ray << " angle "
        //        << is << ": "
        //        << ray2D[0] << "\n" << std::flush;
        //}
    }
}

//-1 on failure
int RunBellhop(std::string FileRoot, void * result) {
    ray2DAll = reinterpret_cast<ray2DPt*>(result); //no checking

    bool singlethread = false;
    if (FileRoot.empty()) {
        std::cout << "Must provide FileRoot as command-line parameter\n";
        return -1;
    }

    FileRoot = std::string("d:/Users/oldst/Documents/joe/pong/temp/") + FileRoot;

    //std::cout << " size " << sizeof(ray2DPt) << "\n" << std::flush;

    setup(FileRoot, PRTFile, RAYFile, ARRFile, Title, fT,
        Bdry, bdinfo, refl, ssp, atten, Pos, Angles, freqinfo, Beam, beaminfo);
    core_setup(PRTFile, fT, Bdry, bdinfo, atten, Angles, freqinfo, Beam);

    rayID = 0;

    std::vector<std::thread> threads;
    uint32_t cores = singlethread ? 1u : std::max(std::thread::hardware_concurrency(), 1u);
    for (uint32_t i = 0; i < cores; ++i) threads.push_back(std::thread(RayWorker));
    for (uint32_t i = 0; i < cores; ++i) threads[i].join();

    ARRFile.close();
    PRTFile.close();

    return rayID;
}

