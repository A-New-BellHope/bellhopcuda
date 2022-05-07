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
 * Topn, Botn: top, bottom normal vector (outward)
 * DistTop, DistBot: distance (normal to bdry) from the ray to top, bottom boundary
*/
template<bool R3D> HOST_DEVICE inline void Distances(const VEC23<R3D> &rayx, 
    const VEC23<R3D> &Topx, const VEC23<R3D> &Botx,
    const VEC23<R3D> &Topn, const VEC23<R3D> &Botn, real &DistTop, real &DistBot)
{
    VEC23<R3D> dTop = rayx - Topx; // vector pointing from top    bdry to ray
    VEC23<R3D> dBot = rayx - Botx; // vector pointing from bottom bdry to ray
    DistTop = -glm::dot(Topn, dTop);
    DistBot = -glm::dot(Botn, dBot);
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
    Origin<false, false> &org, SSPSegState &iSeg, BdryState<false> &bds, BdryType &Bdry,
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
    EvaluateSSP<false, false>(xs, tinit, o, org, ssp, iSeg);
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
    // read by Reflect, is never set in bdinfo and is left alone from
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
template<bool O3D, bool R3D> HOST_DEVICE inline int32_t RayUpdate(
    const rayPt<R3D> &point0, rayPt<R3D> &point1, rayPt<R3D> &point2,
    real &DistEndTop, real &DistEndBot,
    int32_t &iSmallStepCtr, const Origin<O3D, R3D> &org, SSPSegState &iSeg,
    BdryState<O3D> &bds, BdryType &Bdry, const BdryInfo<O3D> *bdinfo, const ReflectionInfo *refl,
    const SSPStructure *ssp, const FreqInfo *freqinfo, const BeamStructure *Beam)
{
    int32_t numRaySteps = 1;
    bool topRefl, botRefl;
    Step<O3D, R3D>(point0, point1, bds, Beam, org, ssp, iSeg, iSmallStepCtr, topRefl, botRefl);
    /*
    if(point0.x == point1.x){
        printf("Ray did not move from (%g,%g), bailing\n", point0.x.x, point0.x.y);
        bail();
    }
    */
    
    VEC23<O3D> x_o = RayToOceanX(point1.x, org);
    VEC23<O3D> t_o = RayToOceanT(point1.t, org);
    GetBdrySeg(x_o, t_o, bds.top, &bdinfo->top, true);
    GetBdrySeg(x_o, t_o, bds.bot, &bdinfo->bot, false);
    
    //TODO
    // Reflections?
    // Tests that ray at step is is inside, and ray at step is+1 is outside
    // to detect only a crossing from inside to outside
    // DistBeg is the distance at point0, which is saved
    // DistEnd is the distance at point1, which needs to be calculated
    Distances(point1.x, bds.top.x, bds.bot.x,
        bds.top.n, bds.bot.n, DistEndTop, DistEndBot);
    
    // LP: Merging these cases is important for GPU performance.
    if(topRefl || botRefl){
        // printf(topRefl ? "Top reflecting\n" : "Bottom reflecting\n");
        const BdryInfoTopBot<O3D> &bdi = topRefl ? bdinfo->top : bdinfo->bot;
        const BdryStateTopBot<O3D> &bdstb = topRefl ? bds.top : bds.bot;
        const HSInfo &hs = topRefl ? Bdry.Top.hs : Bdry.Bot.hs;
        const ReflectionInfoTopBot &refltb = topRefl ? refl->top : refl->bot;
        ReflCurvature<O3D> rcurv;
        VEC23<R3D> tInt(RL(0.0));
        VEC23<O3D> nInt;
        
        if constexpr(O3D){
            // LP: FORTRAN actually checks if the whole string is just "C", not just the first char
            if(bdi.type[0] == 'C'){
                //TODO 2D-3D
                real s1 = (point1.x - bdstb.x) / (bdstb.lSeg.x.max - bdstb.lSeg.x.min);
                real s2 = (point1.y - bdstb.y) / (bdstb.lSeg.y.max - bdstb.lSeg.y.min);
                real m1 = FL(1.0) - s1;
                real m2 = FL(1.0) - s2;
                
                BdryPtFull<true> *bd00 = &bdi.bd[(bdstb.Iseg.x  )*bdi.NPts.y+bdstb.Iseg.y  ];
                BdryPtFull<true> *bd01 = &bdi.bd[(bdstb.Iseg.x  )*bdi.NPts.y+bdstb.Iseg.y+1];
                BdryPtFull<true> *bd10 = &bdi.bd[(bdstb.Iseg.x+1)*bdi.NPts.y+bdstb.Iseg.y  ];
                BdryPtFull<true> *bd11 = &bdi.bd[(bdstb.Iseg.x+1)*bdi.NPts.y+bdstb.Iseg.y+1];
                
                nInt = bd00.Noden * m1 * m2 +
                       bd10.Noden * s1 * m2 +
                       bd11.Noden * s1 * s2 +
                       bd01.Noden * m1 * s2;
                rcurv.z_xx = bd00.z_xx;
                rcurv.z_xy = bd00.z_xy;
                rcurv.z_yy = bd00.z_yy;
                
                rcurv.kappa_xx = bd00.kappa_xx;
                rcurv.kappa_xy = bd00.kappa_xy;
                rcurv.kappa_yy = bd00.kappa_yy;
            }else{
                nInt = bdstb.n;
                rcurv.z_xx = rcurv.z_xy = rcurv.z_yy = FL(0.0);
                rcurv.kappa_xx = rcurv.kappa_xy = rcurv.kappa_yy = FL(0.0);
            }
        }else{
            BdryPtFull<false> *bd0 = &bdi.bd[bdstb.Iseg];
            BdryPtFull<false> *bd1 = &bd0[1]; // LP: next segment
            // LP: FORTRAN actually checks if the whole string is just "C", not just the first char
            if(bdi.type[0] == 'C'){
                real sss = glm::dot(point1.x - bdstb.x, bd0->t) / bd0->Len; // proportional distance along segment
                nInt = (FL(1.0) - sss) * bd0->Noden + sss * bd1->Noden;
                tInt = (FL(1.0) - sss) * bd0->Nodet + sss * bd1->Nodet;
            }else{
                nInt = bd0->n; // normal is constant in a segment
                tInt = bd0->t;
            }
            rcurv.kappa = bd0->kappa;
        }
        
        Reflect<O3D, R3D>(point1, point2, hs, topRefl, tInt, nInt, rcurv,
            freqinfo->freq0, refltb, Beam, org, ssp, iSeg);
        //Incrementing bounce count moved to Reflect
        numRaySteps = 2;
        Distances(point2.x, bds.top.x, bds.bot.x,
            bds.top.n, bds.bot.n, DistEndTop, DistEndBot);
    }
    
    return numRaySteps;
}

/**
 * Has the ray left the box, lost its energy, escaped the boundaries, or 
 * exceeded storage limit?
 * [LP: 2D-3D only:]
 * this should be modified to have a single box
 * no need to test point.x.x, for instance, against several limits; calculate one limit in advance
 * LP: Also updates DistBegTop, DistBegBot.
 */
template<bool O3D, bool R3D> HOST_DEVICE inline bool RayTerminate(
    const rayPt<O3D> &opoint, const rayPt<R3D> &rpoint,
    int32_t &Nsteps, int32_t is, const int32_t &iSmallStepCtr,
    real &DistBegTop, real &DistBegBot, const real &DistEndTop, const real &DistEndBot,
    const Origin<O3D, R3D> &org, const BdryInfo<O3D> *bdinfo, const BeamStructure *Beam
    )
{
    bool leftbox, escapedboundaries, backward, toomanysmallsteps;
    bool escaped0bdry, escapedNbdry;
    if constexpr(O3D){
        leftbox = STD::abs(opoint.x.x - org.xs.x) > Beam->Box.x ||
                  STD::abs(opoint.x.y - org.xs.y) > Beam->Box.y ||
                  STD::abs(opoint.x.z - org.xs.z) > Beam->Box.z;
        escaped0bdry = 
            x.x < bhc::max(bdinfo->bot.bd[0].x.x, bdinfo->top.bd[0].x.x) ||
            x.y < bhc::max(bdinfo->bot.bd[0].x.y, bdinfo->top.bd[0].x.y);
        escapedNbdry =
            x.x > bhc::max(bdinfo->bot.bd[(bdinfo->bot.NPts.x-1)*bdinfo->bot.NPts.y].x.x, 
                           bdinfo->top.bd[(bdinfo->top.NPts.x-1)*bdinfo->top.NPts.y].x.x) ||
            x.y > bhc::max(bdinfo->bot.bd[bdinfo->bot.NPts.y-1].x.y, 
                           bdinfo->top.bd[bdinfo->top.NPts.y-1].x.y);
        escapedboundaries = escaped0bdry || escapedNbdry;
        toomanysmallsteps = iSmallStepCtr > 50;
    }else{
        leftbox = STD::abs(rpoint.x.x) > Beam->Box.r ||
                  STD::abs(rpoint.x.y) > Beam->Box.z;
        escaped0bdry = escapedNbdry = false;
        escapedboundaries = (DistBegTop < FL(0.0) && DistEndTop < FL(0.0)) ||
                            (DistBegBot < FL(0.0) && DistEndBot < FL(0.0));
        toomanysmallsteps = false; // LP: Simply never checked in 2D.
    }
    bool lostenergy = rpoint.Amp < FL(0.005);
    if constexpr(O3D && !R3D){
        backward = rpoint.t.x < FL(0.0); // kills off a backward traveling ray
    }else{
        backward = false; // LP: Commented out for 2D, absent for 3D as would not make sense.
    }
    if(leftbox || lostenergy || escapedboundaries || backward){
        /*
        if(leftbox){
            printf("Ray left beam box (%g,%g)\n", Beam->Box.r, Beam->Box.z);
        }else if(escapedboundaries){
            printf("Ray escaped boundaries DistBegTop %g DistEndTop %g DistBegBot %g DistEndBot %g\n",
                DistBegTop, DistEndTop, DistBegBot, DistEndBot);
        }else if(lostenergy){
            printf("Ray energy dropped to %g\n", point.Amp);
        }else if(backward){
            printf("Ray is going backwards\n");
        }
        */
        if(O3D && escaped0bdry){
            // LP: 2D-3D and 3D: If escapes the boundary only to the negative
            // side, stop without including the current step.
            Nsteps = is;
        }else{
            Nsteps = is + 1;
        }
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
