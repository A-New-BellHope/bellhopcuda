#pragma once
#include "trace.hpp"
#include "influence.hpp"

HOST_DEVICE inline int32_t GetNumJobs(const Position *Pos, const AnglesStructure *Angles)
{
    return Pos->NSz * (Angles->iSingle_alpha >= 1 ? 1 : Angles->Nalpha);
}

/**
 * Returns whether the job should continue.
 * `is` changed to `isrc` because `is` is used for steps
 */
HOST_DEVICE inline bool GetJobIndices(int32_t &isrc, int32_t &ialpha, int32_t job,
    const Position *Pos, const AnglesStructure *Angles)
{
    if(Angles->iSingle_alpha >= 1){
        isrc = job;
        ialpha = Angles->iSingle_alpha - 1; //iSingle_alpha is 1-indexed because how defined in env file
    }else{
        isrc = job / Angles->Nalpha;
        ialpha = job % Angles->Nalpha;
    }
    return (isrc < Pos->NSz);
}

/**
 * Main ray tracing function for ray path output mode.
 */
HOST_DEVICE inline void MainRayMode(int32_t isrc, int32_t ialpha, real &SrcDeclAngle,
    ray2DPt *ray2D, int32_t &Nsteps,
    const BdryType *ConstBdry, const BdryInfo *bdinfo, const ReflectionInfo *refl,
    const SSPStructure *ssp, const Position *Pos, const AnglesStructure *Angles,
    const FreqInfo *freqinfo, const BeamStructure *Beam, const BeamInfo *beaminfo)
{
    real DistBegTop, DistEndTop, DistBegBot, DistEndBot;
    int32_t IsegTop, IsegBot, iSegz, iSegr;
    vec2 gradc, rTopSeg, rBotSeg;
    BdryType Bdry;
    
    if(!RayInit(isrc, ialpha, SrcDeclAngle, ray2D[0], gradc,
        DistBegTop, DistBegBot, IsegTop, IsegBot, rTopSeg, rBotSeg, iSegz, iSegr,
        Bdry, ConstBdry, bdinfo, refl, ssp, Pos, Angles, freqinfo, Beam, beaminfo))
    {
        Nsteps = 1;
        return;
    }
    
    int32_t iSmallStepCtr = 0;
    int32_t is = 0; // index for a step along the ray
    
    for(int32_t istep = 0; istep<MaxN-1; ++istep){
        is += RayUpdate(ray2D[is], ray2D[is+1], ray2D[is+2], 
            DistBegTop, DistBegBot, DistEndTop, DistEndBot,
            IsegTop, IsegBot, rTopSeg, rBotSeg, iSmallStepCtr, iSegz, iSegr,
            Bdry, bdinfo, refl, ssp, freqinfo, Beam);
        if(RayTerminate(ray2D[is], Nsteps, is, DistBegTop, DistBegBot,
            DistEndTop, DistEndBot, Beam)) break;
        if(Nsteps >= 0 && is > Nsteps){
            Nsteps = is + 1;
            break;
        }
    }
}

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

/**
 * LP: Write TL results
 */
inline void FinalizeTLMode(cpxf *&uAllSources, DirectOFile &SHDFile,
    const SSPStructure *ssp, const Position *Pos, const AnglesStructure *Angles,
    const FreqInfo *freqinfo, const BeamStructure *Beam)
{
    for(int32_t isrc=0; isrc<Pos->NSz; ++isrc){
        cpx ccpx;
        int32_t iSegz = 0, iSegr = 0;
        EvaluateSSPCOnly(vec2(RL(0.0), Pos->Sz[isrc]), vec2(RL(1.0), RL(0.0)),
            ccpx, freqinfo->freq0, ssp, iSegz, iSegr);
        ScalePressure(Angles->Dalpha, ccpx.real(), Pos->Rr, 
            &uAllSources[isrc * Pos->NRz_per_range * Pos->NRr], 
            Pos->NRz_per_range, Pos->NRr, Beam->RunType, freqinfo->freq0);
        int32_t IRec = 10 + Pos->NRz_per_range * isrc;
        for(int32_t Irz1 = 0; Irz1 < Pos->NRz_per_range; ++Irz1){
            SHDFile.rec(IRec);
            for(int32_t r=0; r < Pos->NRr; ++r){
                DOFWRITEV(SHDFile, uAllSources[(isrc * Pos->NRz_per_range + Irz1) * Pos->NRr + r]);
            }
            ++IRec;
        }
    }
    deallocate(uAllSources);
}

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
