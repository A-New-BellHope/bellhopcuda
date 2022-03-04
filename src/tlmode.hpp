#pragma once
#include "common.hpp"
#include "trace.hpp"
#include "ssp.hpp"
#include "influence.hpp"

/**
 * for a TL calculation, allocate space for the pressure matrix
 */
inline void InitTLMode(cpxf *&uAllSources, const Position *Pos,
    const BeamStructure *Beam)
{
    size_t n = Pos->NSz * Pos->NRz_per_range * Pos->NRr;
    uAllSources = allocate<cpxf>(n);
    memset(uAllSources, 0, n * sizeof(cpxf));
}

void FinalizeTLMode(std::string FileRoot, const bhcParams &params, bhcOutputs &outputs);

/**
 * Main ray tracing function for TL, eigen, and arrivals runs.
 */
HOST_DEVICE inline void MainFieldModes(int32_t isrc, int32_t ialpha, real &SrcDeclAngle,
    cpxf *uAllSources,
    const BdryType *ConstBdry, const BdryInfo *bdinfo, const ReflectionInfo *refl,
    const SSPStructure *ssp, const Position *Pos, const AnglesStructure *Angles,
    const FreqInfo *freqinfo, const BeamStructure *Beam, const BeamInfo *beaminfo,
    EigenInfo *eigen, const ArrInfo *arrinfo)
{
    real DistBegTop, DistEndTop, DistBegBot, DistEndBot;
    int32_t IsegTop, IsegBot, iSegz, iSegr;
    vec2 gradc, rTopSeg, rBotSeg;
    BdryType Bdry;
    
    ray2DPt point0, point1, point2;
    InfluenceRayInfo inflray;
    
    if(!RayInit(isrc, ialpha, SrcDeclAngle, point0, gradc,
        DistBegTop, DistBegBot, IsegTop, IsegBot, rTopSeg, rBotSeg, iSegz, iSegr,
        Bdry, ConstBdry, bdinfo, refl, ssp, Pos, Angles, freqinfo, Beam, beaminfo)) return;
    
    Init_Influence(inflray, point0, isrc, ialpha, Angles->alpha[ialpha], gradc,
        Pos, ssp, iSegz, iSegr, Angles, freqinfo, Beam);
    
    cpxf *u = &uAllSources[isrc * Pos->NRz_per_range * Pos->NRr];
    int32_t iSmallStepCtr = 0;
    int32_t is = 0; // index for a step along the ray
    int32_t Nsteps = 0; // not actually needed in TL mode, debugging only
    
    for(int32_t istep = 0; istep<MaxN-1; ++istep){
        int32_t dStep = RayUpdate(point0, point1, point2, 
            DistBegTop, DistBegBot, DistEndTop, DistEndBot,
            IsegTop, IsegBot, rTopSeg, rBotSeg, iSmallStepCtr, iSegz, iSegr,
            Bdry, bdinfo, refl, ssp, freqinfo, Beam);
        if(!Step_Influence(point0, point1, inflray, is, u, 
            ConstBdry, ssp, iSegz, iSegr, Pos, Beam, eigen, arrinfo)){
            //printf("Step_Influence terminated ray\n");
            break;
        }
        ++is;
        if(dStep == 2){
            if(!Step_Influence(point1, point2, inflray, is, u, 
                ConstBdry, ssp, iSegz, iSegr, Pos, Beam, eigen, arrinfo)) break;
            point0 = point2;
            ++is;
        }else if(dStep == 1){
            point0 = point1;
        }else{
            printf("Invalid dStep: %d\n", dStep);
            bail();
        }
        if(RayTerminate(point0, Nsteps, is, DistBegTop, DistBegBot,
            DistEndTop, DistEndBot, Beam)) break;
    }
    
    //printf("Nsteps %d\n", Nsteps);
}
