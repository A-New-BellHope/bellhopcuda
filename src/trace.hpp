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
 * Copy only specific data from the HSInfo (halfspace info) struct.
 * [LP: FORTRAN] compiler is not accepting the copy of the whole structure at once ...
 * LP: maybe this means actually the whole struct should be copied, but he
 * only copied the elements which were needed?
 */
HOST_DEVICE inline void CopyHSInfo(HSInfo &b, const HSInfo &a)
{
    b.cp  = a.cp;
    b.cs  = a.cs;
    b.rho = a.rho;
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
    if(atiType[1] == 'L') CopyHSInfo(Bdry.Top.hs, Top[IsegTop].hs);
    if(btyType[1] == 'L') CopyHSInfo(Bdry.Bot.hs, Bot[IsegBot].hs);
    
    // Trace the beam (note that Reflect alters the step index is)
    is = 0;
    Distances2D(ray2D[0].x, Top[IsegTop].x, Bot[IsegBot].x, dEndTop, dEndBot,
        Top[IsegTop].n, Bot[IsegBot].n, DistBegTop, DistBegBot);
    
    if(DistBegTop <= 0 || DistBegBot <= 0){
        Beam->Nsteps = 1;
        printf("Terminating the ray trace because the source is on or outside the boundaries\n");
        return; // source must be within the medium
    }
    
    for(int32_t istep = 0; istep<MaxN-1; ++istep){
        is  = istep + 1;
        is1 = istep + 1;
        
        Step2D(ray2D[is], &ray2D[is1], Top[IsegTop].x, Top[IsegTop].n,
            Bot[IsegBot].x, Bot[IsegBot].n, freq, Beam, ssp, iSegz, iSegr, iSmallStepCtr);
        
        // New altimetry segment?
        if(ray2D[is1].x.x < rTopSeg.x || ray2D[is1].x.x > rTopSeg.y){
            GetTopSeg(ray2D[is1].x.x, IsegTop, rTopseg);
            if(atiType[1] == 'L') CopyHSInfo(Bdry.Top.hs, Top[IsegTop].hs); // grab the geoacoustic info for the new segment
        }
        
        // New bathymetry segment?
        if(ray2D[is1].x.x < rBotSeg.x || ray2D[is1].x.x > rBotSeg.y){
            GetBotSeg(ray2D[is1].x.x, IsegBot, rBotseg);
            if(btyType[1] == 'L') CopyHSInfo(Bdry.Bot.hs, Bot[IsegBot].hs); // grab the geoacoustic info for the new segment
        }
        
        // Reflections?
        // Tests that ray at step is is inside, and ray at step is+1 is outside
        // to detect only a crossing from inside to outside
        // DistBeg is the distance at step is,   which is saved
        // DistEnd is the distance at step is+1, which needs to be calculated
        
        Distances2D(ray2D[is1].x, Top[IsegTop].x, Bot[IsegBot].x, dEndTop, dEndBot,
            Top[IsegTop].n, Bot[IsegBot].n, DistEndTop, DistEndBot);
        
        if(DistBegTop > RC(0.0) && DistEndTop <= RC(0.0)){ // test top reflection
            
            if(atiType[0] == 'C'){ // LP: Actually checking if the whole string is just "C", not just the first char
                sss = glm::dot(dEndTop, Top[IsegTop].t) / Top[IsegTop].Len; // proportional distance along segment
                TopnInt = (RC(1.0) - sss) * Top[IsegTop].Noden + sss * Top[IsegTop+1].Noden;
                ToptInt = (RC(1.0) - sss) * Top[IsegTop].Nodet + sss * Top[IsegTop+1].Nodet;
            }else{
                TopnInt = Top[IsegTop].n; // normal is constant in a segment
                ToptInt = Top[IsegTop].t;
            }
            
            Reflect2D(is, Bdry.Top.hs, true, ToptInt, TopnInt, Top[IsegTop].Kappa, 
                RTop, NTopPTS);
            ++ray2D[is+1].NumTopBnc;
            
            Distances2D(ray2D[is+1].x, Top[IsegTop].x, Bot[IsegBot].x, dEndTop, dEndBot,
                Top[IsegTop].n, Bot[IsegBot].n, DistEndTop, DistEndBot);
                
        }else if(DistBegBot > RC(0.0) && DistEndBot <= RC(0.0)){ // test bottom reflection
            
            if(btyType[0] == 'C'){ // LP: Actually checking if the whole string is just "C", not just the first char
                sss = glm::dot(dEndBot, Bot[IsegBot].t) / Bot[IsegBot].Len; // proportional distance along segment
                BotnInt = (RC(1.0) - sss) * Bot[IsegBot].Noden + sss * Bot[IsegBot+1].Noden;
                BottInt = (RC(1.0) - sss) * Bot[IsegBot].Nodet + sss * Bot[IsegBot+1].Nodet;
            }else{
                BotnInt = Bot[IsegBot].n; // normal is constant in a segment
                BottInt = Bot[IsegBot].t;
            }
            
            Reflect2D(is, Bdry.Bot.hs, true, BottInt, BotnInt, Bot[IsegBot].Kappa, 
                RBot, NBotPTS);
            ++ray2D[is+1].NumBotBnc;
            
            Distances2D(ray2D[is+1].x, Top[IsegTop].x, Bot[IsegBot].x, dEndTop, dEndBot,
                Top[IsegTop].n, Bot[IsegBot].n, DistEndTop, DistEndBot);
            
        }
        
        // Has the ray left the box, lost its energy, escaped the boundaries, or exceeded storage limit?
        if( STD::abs(ray2D[is+1].x.x) > Beam->Box.r ||
            STD::abs(ray2d[is+1].x.y) > Beam->Box.z ||
            ray2d[is+1].Amp < RC(0.005) ||
            (DistBegTop < RC(0.0) && DistEndTop < RC(0.0)) ||
            (DistBegBot < RC(0.0) && DistEndBot < RC(0.0)) 
            // ray2d[is+1].t.x < 0 // this last test kills off a backward traveling ray
        ){
            Beam->Nsteps = is + 1;
            break;
        }else if(is >= MaxN - 4){
            printf("Warning in TraceRay2D: Insufficient storage for ray trajectory\n");
            Beam->Nsteps = is;
            break;
        }
        
        DistBegTop = DistEndTop;
        DistBegBot = DistEndBot;        
    }
}
