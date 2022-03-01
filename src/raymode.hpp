#pragma once
#include "common.hpp"

void OpenRAYFile(LDOFile &RAYFile, std::string FileRoot, bool ThreeD, 
    const bhcParams &params);

void WriteRay2D(real alpha0, int32_t Nsteps1, LDOFile &RAYFile,
    const BdryType *Bdry, ray2DPt *ray2D);
    
void InitRayMode(RayInfo *rayinfo, const bhcParams &params);
void FinalizeRayMode(RayInfo *rayinfo, std::string FileRoot, const bhcParams &params);

inline bool IsRayCopyMode(const RayInfo *rayinfo)
{
    return (size_t)rayinfo->MaxPoints < (size_t)MaxN * (size_t)rayinfo->NRays;
}

inline bool RunRay(RayInfo *rayinfo, const bhcParams &params, ray2DPt *localmem,
    int32_t job, int32_t isrc, int32_t ialpha)
{
    ray2DPt *ray2D;
    if(IsRayCopyMode(rayinfo)){
        ray2D = localmem;
    }else{
        ray2D = &rayinfo->raymem[job * MaxN];
    }
    memset(ray2D, 0xFE, MaxN * sizeof(ray2DPt)); //Set to garbage values for debugging
    
    real SrcDeclAngle;
    int32_t Nsteps = -1;
    MainRayMode(isrc, ialpha, SrcDeclAngle, ray2D, Nsteps,
        params.Bdry, params.bdinfo, params.refl, params.ssp, params.Pos,
        params.Angles, params.freqinfo, params.Beam, params.beaminfo);
    
    bool ret = true;
    if(IsRayCopyMode(rayinfo)){
        uint32_t p = AtomicFetchAdd(&rayinfo->NPoints, Nsteps);
        if(p + Nsteps > rayinfo->MaxPoints){
            std::cout << "Ran out of memory for rays\n";
            rayinfo->results[job].ray2D = nullptr;
            ret = false;
        }else{
            rayinfo->results[job].ray2D = &rayinfo->raymem[p];
            memcpy(rayinfo->results[job].ray2D, localmem, Nsteps * sizeof(ray2DPt));
        }
    }else{
        rayinfo->results[job].ray2D = ray2D;
    }
    rayinfo->results[job].SrcDeclAngle = SrcDeclAngle;
    rayinfo->results[job].Nsteps = Nsteps;
    
    return ret;
}
