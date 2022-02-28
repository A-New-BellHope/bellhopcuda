//Joe Snider
//12/21
//
//Testing shared library target for the bellhop translation.
//Not real general.

#define BELLHOPCXX_EXPORTS

#include "trace_lib.hpp"
#include "setup.hpp"
#include "main.hpp"
#include "output.hpp"

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
DirectOFile SHDFile;
EigenInfo* eigen;
ArrInfo* arrinfo;

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

/// <summary>
/// start on threads for a ray calculation. Mostly copied from tl_cxx.
/// </summary>
void RayWorker()
{
    while (true) {
        int32_t ray;
        //ray = rayID++; //should work, but something is strange with operator++ and atomics
        //not much cost to mutex the job number (ray).
        {
            std::lock_guard<std::mutex> g(rayfilemutex);
            ray = rayID;
            rayID += 1;
        }
        ray2DPt* ray2D = ray2DAll + MaxN * ray;
        int32_t isrc, ialpha;
        if (!GetJobIndices(isrc, ialpha, ray, Pos, Angles)) break;

        memset(ray2D, 0xFE, MaxN * sizeof(ray2DPt)); //Set to garbage values for debugging

        real SrcDeclAngle;
        int32_t Nsteps = -1;
        MainRayMode(isrc, ialpha, SrcDeclAngle, ray2D, Nsteps,
            Bdry, bdinfo, refl, ssp, Pos, Angles, freqinfo, Beam, beaminfo);
    }
}

int RunBellhop(const char* cFileRoot, void * result) {
    ray2DAll = reinterpret_cast<ray2DPt*>(result); //no checking

    bool singlethread = false;
    std::string FileRoot = cFileRoot;
    if (FileRoot.empty()) {
        std::cout << "Must provide FileRoot as command-line parameter\n";
        return -1;
    }

    //FileRoot = std::string("d:/Users/oldst/Documents/joe/pong/temp/") + FileRoot;

    setup(FileRoot, PRTFile, RAYFile, SHDFile, Title, fT,
        Bdry, bdinfo, refl, ssp, atten, Pos, Angles, freqinfo, Beam, beaminfo, eigen, arrinfo);
    InitCommon(Pos, Beam);

    rayID = 0;

    std::vector<std::thread> threads;
    uint32_t cores = singlethread ? 1u : std::max(std::thread::hardware_concurrency(), 1u);
    for (uint32_t i = 0; i < cores; ++i) threads.push_back(std::thread(RayWorker));
    for (uint32_t i = 0; i < cores; ++i) threads[i].join();

    ARRFile.close();
    PRTFile.close();

    return rayID;
}
