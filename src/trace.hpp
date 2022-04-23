/*
bellhopcxx / bellhopcuda - C++/CUDA port of BELLHOP underwater acoustics simulator
Copyright (C) 2021-2022 The Regents of the University of California
c/o Jules Jaffe team at SIO / UCSD, jjaffe@ucsd.edu
Based on BELLHOP, which is Copyright (C) 1983-2020 Michael B. Porter

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <https://www.gnu.org/licenses/>.
*/
#pragma once
#include "common.hpp"
#include "step.hpp"
#include "reflect.hpp"

namespace bhc {

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
    b.cP  = a.cP;
    b.cS  = a.cS;
    b.rho = a.rho;
}

/**
 * LP: Pulled out ray update loop initialization. Returns whether to continue
 * with the ray trace. Only call for valid ialpha w.r.t. Angles->iSingleAlpha.
 * Original comments follow.
 * 
 * DistBegTop etc.: Distances from ray beginning, end to top and bottom
 */
HOST_DEVICE inline bool RayInit(int32_t isrc, int32_t ialpha, real &SrcDeclAngle,
    rayPt<false> &point0, vec2 &gradc, real &DistBegTop, real &DistBegBot, 
    SSPSegState &iSeg, BdryState<false> &bds, BdryType &Bdry,
    const BdryType *ConstBdry, const BdryInfo<false> *bdinfo,
    const SSPStructure *ssp, const Position *Pos, const AnglesStructure *Angles,
    const FreqInfo *freqinfo, const BeamStructure *Beam, const BeamInfo *beaminfo)
{
    if(isrc < 0 || ialpha < 0 || isrc >= Pos->NSz || ialpha >= Angles->Nalpha){
        printf("Invalid isrc %d ialpha %d\n", isrc, ialpha);
        bail();
    }
    
    // LP: This part from BellhopCore
    float omega = FL(2.0) * REAL_PI * freqinfo->freq0;
    vec2 xs = vec2(FL(0.0), Pos->Sz[isrc]); // x-y coordinate of the source
    
    // LP: Changed in BELLHOP to just reinitialize before every ray.
    iSeg.z = iSeg.r = 0;
    
    real alpha = Angles->alpha[ialpha]; // initial angle
    SrcDeclAngle = RadDeg * alpha; // take-off declination angle in degrees
    
    SSPOutputs<false> o;
    vec2 tinit = vec2(STD::cos(alpha), STD::sin(alpha));
    EvaluateSSP<false, false>(xs, tinit, o, ssp, iSeg);
    gradc = o.gradc;
    
    // Are there enough beams?
    real DalphaOpt = STD::sqrt(o.ccpx.real() / (FL(6.0) * freqinfo->freq0 * Pos->Rr[Pos->NRr-1]));
    int32_t NalphaOpt = 2 + (int)((Angles->alpha[Angles->Nalpha-1] - Angles->alpha[0]) / DalphaOpt);
    
    if(Beam->RunType[0] == 'C' && Angles->Nalpha < NalphaOpt && ialpha == 0){
        printf("Warning in " BHC_PROGRAMNAME " : Too few beams\nNalpha should be at least = %d\n", NalphaOpt);
    }
    
    int32_t ibp = BinarySearchLEQ(beaminfo->SrcBmPat, beaminfo->NSBPPts, 2, 0, SrcDeclAngle);
    ibp = bhc::min(ibp, beaminfo->NSBPPts-2); // don't go past end of table
    
    // linear interpolation to get amplitude
    real s = (SrcDeclAngle - beaminfo->SrcBmPat[2*ibp+0]) / (beaminfo->SrcBmPat[2*(ibp+1)+0] - beaminfo->SrcBmPat[2*ibp+0]);
    float Amp0 = (FL(1.0) - s) * beaminfo->SrcBmPat[2*ibp+1] + s * beaminfo->SrcBmPat[2*(ibp+1)+1]; // initial amplitude
    
    // Lloyd mirror pattern for semi-coherent option
    if(Beam->RunType[0] == 'S')
        Amp0 *= STD::sqrt(FL(2.0)) * STD::abs(STD::sin(omega / o.ccpx.real() * xs.y * STD::sin(alpha)));
        
    // LP: This part from TraceRay<false>
    point0.NumTopBnc = 0;
    point0.NumBotBnc = 0;
    point0.x         = xs;
    point0.t         = tinit / o.ccpx.real();
    point0.p         = vec2(FL(1.0), FL(0.0));
    point0.q         = vec2(FL(0.0), FL(1.0));
    point0.c         = o.ccpx.real();
    point0.Amp       = Amp0;
    point0.Phase     = FL(0.0);
    point0.tau       = cpx(FL(0.0), FL(0.0));
    
    // second component of qv is not used in geometric beam tracing
    // set I.C. to 0 in hopes of saving run time
    if(Beam->RunType[1] == 'G') point0.q = vec2(FL(0.0), FL(0.0));
    
    bds.top.Iseg = 0; bds.bot.Iseg = 0;
    GetBdrySeg(xs, point0.t, bds.top, &bdinfo->top, true ); // identify the top    segment above the source
    GetBdrySeg(xs, point0.t, bds.bot, &bdinfo->bot, false); // identify the bottom segment below the source
    
    // convert range-dependent geoacoustic parameters from user to program units
    // LP: BELLHOP uses all values from ConstBdry except replaces cP, cS, and
    // rho from the current segment in bdinfo. rho is read from files in
    // ReadBoundary, and cP and cS are computed in core_setup. bc, which is also
    // read by Reflect2D, is never set in bdinfo and is left alone from
    // ConstBdry.
    Bdry = *ConstBdry;
    if(bdinfo->top.type[1] == 'L') CopyHSInfo(Bdry.Top.hs, bdinfo->top.bd[bds.top.Iseg].hs);
    if(bdinfo->bot.type[1] == 'L') CopyHSInfo(Bdry.Bot.hs, bdinfo->bot.bd[bds.bot.Iseg].hs);
    // printf("btyType cP top bot %c%c (%g,%g) (%g,%g)\n", bdinfo->bot.type[0], bdinfo->bot.type[1],
    //     Bdry.Top.hs.cP.real(), Bdry.Top.hs.cP.imag(),
    //     Bdry.Bot.hs.cP.real(), Bdry.Bot.hs.cP.imag());
    
    vec2 dEndTop_temp, dEndBot_temp;
    Distances2D(point0.x, bds.top.x, bds.bot.x, dEndTop_temp, dEndBot_temp, 
        bds.top.n, bds.bot.n, DistBegTop, DistBegBot);
        
    if(DistBegTop <= FL(0.0) || DistBegBot <= FL(0.0)){
        printf("Terminating the ray trace because the source is on or outside the boundaries\n");
        // printf("xs (%g,%g) Bot.x (%g,%g) Bot.n (%g,%g) DistBegBot %g\n",
        //     xs.x, xs.y, bds.bot.x.x, bds.bot.x.y,
        //     bds.bot.n.x, bds.bot.n.y, DistBegBot);
        return false; // source must be within the medium
    }
    
    return true;
}

/**
 * LP: Pulled out contents of ray update loop. Returns the number of ray steps
 * taken (i.e. normally 1, or 2 if reflected).
 */
HOST_DEVICE inline int32_t RayUpdate(
    const rayPt<false> &point0, rayPt<false> &point1, rayPt<false> &point2,
    real &DistEndTop, real &DistEndBot,
    int32_t &iSmallStepCtr, SSPSegState &iSeg,
    BdryState<false> &bds, BdryType &Bdry, const BdryInfo<false> *bdinfo, const ReflectionInfo *refl,
    const SSPStructure *ssp, const FreqInfo *freqinfo, const BeamStructure *Beam)
{
    int32_t numRaySteps = 1;
    bool topRefl, botRefl;
    Step<false>(point0, point1, bds, Beam, ssp, iSeg, iSmallStepCtr, topRefl, botRefl);
    /*
    if(point0.x == point1.x){
        printf("Ray did not move from (%g,%g), bailing\n", point0.x.x, point0.x.y);
        bail();
    }
    */
    
    // New altimetry segment?
    if(    point1.x.x < bds.top.lSeg.min || (point1.x.x == bds.top.lSeg.min && point1.t.x < FL(0.0))
        || point1.x.x > bds.top.lSeg.max || (point1.x.x == bds.top.lSeg.max && point1.t.x >= FL(0.0))){
        GetBdrySeg(point1.x, point1.t, bds.top, &bdinfo->top, true);
        if(bdinfo->top.type[1] == 'L') CopyHSInfo(Bdry.Top.hs, bdinfo->top.bd[bds.top.Iseg].hs); // grab the geoacoustic info for the new segment
    }
    
    // New bathymetry segment?
    if(    point1.x.x < bds.bot.lSeg.min || (point1.x.x == bds.bot.lSeg.min && point1.t.x < FL(0.0))
        || point1.x.x > bds.bot.lSeg.max || (point1.x.x == bds.bot.lSeg.max && point1.t.x >= FL(0.0))){
        GetBdrySeg(point1.x, point1.t, bds.bot, &bdinfo->bot, false);
        if(bdinfo->bot.type[1] == 'L') CopyHSInfo(Bdry.Bot.hs, bdinfo->bot.bd[bds.bot.Iseg].hs); // grab the geoacoustic info for the new segment
    }
    
    // Reflections?
    // Tests that ray at step is is inside, and ray at step is+1 is outside
    // to detect only a crossing from inside to outside
    // DistBeg is the distance at point0, which is saved
    // DistEnd is the distance at point1, which needs to be calculated
    vec2 dEndTop, dEndBot;
    Distances2D(point1.x, bds.top.x, bds.bot.x, dEndTop, dEndBot,
        bds.top.n, bds.bot.n, DistEndTop, DistEndBot);
    
    // LP: Merging these cases is important for GPU performance.
    if(topRefl || botRefl){
        // printf(topRefl ? "Top reflecting\n" : "Bottom reflecting\n");
        const BdryInfoTopBot<false> &bdi = topRefl ? bdinfo->top : bdinfo->bot;
        const BdryStateTopBot<false> &bdstb = topRefl ? bds.top : bds.bot;
        vec2 dEnd = topRefl ? dEndTop : dEndBot;
        BdryPtFull<false> *bd0 = &bdi.bd[bdstb.Iseg];
        BdryPtFull<false> *bd1 = &bd0[1]; // LP: next segment
        vec2 nInt, tInt;
        // LP: FORTRAN actually checks if the whole string is just "C", not just the first char
        if(bdi.type[0] == 'C'){
            real sss = glm::dot(dEnd, bd0->t) / bd0->Len; // proportional distance along segment
            nInt = (FL(1.0) - sss) * bd0->Noden + sss * bd1->Noden;
            tInt = (FL(1.0) - sss) * bd0->Nodet + sss * bd1->Nodet;
        }else{
            nInt = bd0->n; // normal is constant in a segment
            tInt = bd0->t;
        }
        Reflect2D(point1, point2, 
            topRefl ? Bdry.Top.hs : Bdry.Bot.hs,
            topRefl, tInt, nInt, bd0->kappa, freqinfo->freq0,
            topRefl ? refl->RTop : refl->RBot,
            topRefl ? refl->NTopPts : refl->NBotPts,
            Beam, ssp, iSeg);
        //Incrementing bounce count moved to Reflect2D
        numRaySteps = 2;
        Distances2D(point2.x, bds.top.x, bds.bot.x, dEndTop, dEndBot,
            bds.top.n, bds.bot.n, DistEndTop, DistEndBot);
    }
    
    return numRaySteps;
}

/**
 * Has the ray left the box, lost its energy, escaped the boundaries, or 
 * exceeded storage limit?
 * LP: Also updates DistBegTop, DistBegBot.
 */
HOST_DEVICE inline bool RayTerminate(const rayPt<false> &point, int32_t &Nsteps, int32_t is,
    real &DistBegTop, real &DistBegBot, const real &DistEndTop, const real &DistEndBot,
    const BeamStructure *Beam
    )
{
    bool leftbox = STD::abs(point.x.x) > Beam->Box.r ||
                   STD::abs(point.x.y) > Beam->Box.z;
    bool lostenergy = point.Amp < FL(0.005);
    bool escapedboundaries = (DistBegTop < FL(0.0) && DistEndTop < FL(0.0)) ||
                             (DistBegBot < FL(0.0) && DistEndBot < FL(0.0));
    //bool backward = point.t.x < 0; // this last test kills off a backward traveling ray
    if(leftbox || lostenergy || escapedboundaries // || backward
    ){
        /*
        if(leftbox){
            printf("Ray left beam box (%g,%g)\n", Beam->Box.r, Beam->Box.z);
        }else if(lostenergy){
            printf("Ray energy dropped to %g\n", point.Amp);
        }else{
            printf("Ray escaped boundaries DistBegTop %g DistEndTop %g DistBegBot %g DistEndBot %g\n",
                DistBegTop, DistEndTop, DistBegBot, DistEndBot);
        }
        */
        Nsteps = is + 1;
        return true;
    }else if(is >= MaxN - 3){
        printf("Warning in TraceRay: Insufficient storage for ray trajectory\n");
        Nsteps = is;
        return true;
    }
    
    DistBegTop = DistEndTop;
    DistBegBot = DistEndBot;
    return false;
}

}
