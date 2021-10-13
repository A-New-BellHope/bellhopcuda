#pragma once
#include "step.hpp"
#include "boundaries.hpp"

/**
 * Calculates the distances to the boundaries
 * Formula differs from JKPS because code uses outward pointing normals
 * 
 * rayx: ray coordinate
 * Topx, Botx: top, bottom coordinate
 * dTop, dBot: vector pointing from top, bottom bdry to ray
 * Topn, Botn: top, bottom normal vector (outward)
 * DistTop, DistBot: distance (normal to bdry) from the ray to top, bottom boundary
*/
HOST_DEVICE inline void Distances2D(const vec2 &rayx, 
    const vec2 &Topx, const vec2 &Botx, vec2 &dTop, vec2 &dBot,
    const vec2 &Topn, const vec2 &Botn, real &DistTop, real &DistBot)
{
    dTop = rayx - Topx; // vector pointing from top    to ray
    dBot = rayx - Botx; // vector pointing from bottom to ray
    DistTop = -glm::dot(Topn, dTop);
    DistBot = -glm::dot(Botn, dBot);
}

/**
 * Traces the beam corresponding to a particular take-off angle
 * 
 * xs: x-y coordinate of the source
 * alpha: initial angle
 * Amp0: initial amplitude
 */
HOST_DEVICE inline void TraceRay2D(vec2 xs, real alpha, real Amp0,
    const real &freq, BeamStructure *Beam, const SSPStructure *ssp,
    const BdryPtFull *Top, const BdryPtFull *Bot,
    const char *atiType, const char *btyType,
    ray2DPt *ray2D)
{
    int32_t is, is1; // index for a step along the ray
    cpx ccpx;
    vec2 gradc, dEndTop, dEndBot, TopnInt, BotnInt, ToptInt, BottInt;
    real crr, crz, czz, rho;
    real DistBegTop, DistEndTop, DistBegBot, DistEndBot; // Distances from ray beginning, end to top and bottom
    real sss;
    
    int32_t IsegTop, IsegBot; // indices that point to the current active segment
    vec2 rTopseg, rBotseg; // range intervals defining the current active segment
    
    BdryType Bdry;
    
    // Initial conditions
    
    int32_t iSmallStepCtr = 0, iSegz = 0, iSegr = 0;
    EvaluateSSP(xs, ccpx, gradc, crr, crz, czz, rho, freq, iSegz, iSegr);
    ray2D[0] = {
        .c         = ccpx.real(),
        .x         = xs,
        .t         = vec2(STD::cos(alpha), STD::sin(alpha)) / ccpx.real(),
        .p         = vec2(RC(1.0), RC(0.0)),
        .q         = vec2(RC(0.0), RC(1.0)),
        .tau       = cpx(RC(0.0), RC(0.0)),
        .Amp       = Amp0,
        .Phase     = RC(0.0),
        .NumTopBnc = 0,
        .NumBotBnc = 0
    };
    
    // second component of qv is not used in geometric beam tracing
    // set I.C. to 0 in hopes of saving run time
    if(Beam->RunType[1] == 'G') ray2D[0].q = vec2(RC(0.0), RC(0.0));
    
    GetTopSeg(xs.x, IsegTop, rTopseg); // identify the top    segment above the source
    GetBotSeg(xs.x, IsegBot, rBotseg); // identify the bottom segment below the source
    
    // convert range-dependent geoacoustic parameters from user to program units
    // [LP: FORTRAN] compiler is not accepting the copy of the whole structure at once ...
    // LP: maybe this means actually the whole struct should be copied, but he
    // only copied the elements which were needed?
    if(atiType[1] == 'L'){
        Bdry.Top.hs.cp  = Top[IsegTop].hs.cp;
        Bdry.Top.hs.cs  = Top[IsegTop].hs.cs;
        Bdry.Top.hs.rho = Top[IsegTop].hs.rho;
    }
    
    if(btyType[1] == 'L'){
        Bdry.Bot.hs.cp  = Bot[IsegBot].hs.cp;
        Bdry.Bot.hs.cs  = Bot[IsegBot].hs.cs;
        Bdry.Bot.hs.rho = Bot[IsegBot].hs.rho;
    }
    
    // Trace the beam (note that Reflect alters the step index is)
    is = 0;
    Distances2D(ray2D[0].x, Top[IsegTop].x, Bot[IsegBot].x, dEndTop, dEndBot,
        Top[IsegTop].n, Bot[IsegBot].n, DistBegTop, DistBegBot);
    
    if(DistBegTop <= 0 || DistBegBot <= 0){
        Beam->Nsteps = 1;
        printf("Terminating the ray trace because the source is on or outside the boundaries\n");
        return; // source must be within the medium
    }
    
    for(int32_t istep = 0; istep<
}
