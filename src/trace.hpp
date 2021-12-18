#pragma once
#include "common.hpp"
#include "boundary.hpp"
#include "refcoef.hpp"
#include "ssp.hpp"
#include "sourcereceiver.hpp"
#include "angles.hpp"
#include "beams.hpp"
#include "step.hpp"

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
 * LP: No function description given.
 * 
 * hs: half-space properties
 * isTop: Flag indicating bottom or top reflection
 * tBdry, nBdry: Tangent and normal to the boundary
 * kappa: Boundary curvature
 * RefC: reflection coefficient
 */
HOST_DEVICE inline void Reflect2D(const ray2DPt &oldPoint, ray2DPt &newPoint,
    const HSInfo &hs, bool isTop,
    const vec2 &tBdry, const vec2 &nBdry, real kappa, real freq,
    const ReflectionCoef *RefC, int32_t Npts,
    const BeamStructure *Beam, const SSPStructure *ssp,
    int32_t &iSegz, int32_t &iSegr)
{
    cpx ccpx;
    vec2 gradc;
    real crr, crz, czz, rho; // derivatives of sound speed
    real rm, rn, Tg, Th;
    vec2 rayt, rayn, rayt_tilde, rayn_tilde;
    real cnjump, csjump; // for curvature change
    real ck, co, si, cco, ssi, pdelta, rddelta, sddelta, theta_bot; // for beam shift
    cpx kx, kz, kzP, kzS, kzP2, kzS2, mu, f, g, y2, y4, Refl; // for tabulated reflection coef.
    cpx ch, a, b, d, sb, delta, ddelta; // for beam shift
    ReflectionCoef RInt;
    real omega = RC(2.0) * M_PI * freq;
    
    Tg = glm::dot(oldPoint.t, tBdry); // component of ray tangent, along boundary
    Th = glm::dot(oldPoint.t, nBdry); // component of ray tangent, normal to boundary
    
    // LP: Incrementing bounce count moved here
    newPoint.NumTopBnc = oldPoint.NumTopBnc + (isTop ? 1 : 0);
    newPoint.NumBotBnc = oldPoint.NumBotBnc + (isTop ? 0 : 1);
    newPoint.x         = oldPoint.x;
    newPoint.t         = oldPoint.t - RC(2.0) * Th * nBdry; // changing the ray direction
    
    // Calculate the change in curvature
    // Based on formulas given by Muller, Geoph. J. R.A.S., 79 (1984).
    
    EvaluateSSP(oldPoint.x, ccpx, gradc, crr, crz, czz, rho, freq, ssp, iSegz, iSegr); // just to get c [LP: ccpx.real()]
    
    // incident unit ray tangent and normal
    rayt = ccpx.real() * oldPoint.t; // unit tangent to ray
    rayn = vec2(-rayt.y, rayt.x);     // unit normal  to ray
    
    // reflected unit ray tangent and normal (the reflected tangent, normal system has a different orientation)
    rayt_tilde = ccpx.real() * newPoint.t;         // unit tangent to ray
    rayn_tilde = -vec2(-rayt_tilde.y, rayt_tilde.x); // unit normal  to ray
    
    rn = RC(2.0) * kappa / SQ(ccpx.real()) / Th; // boundary curvature correction
    
    // get the jumps (this could be simplified, e.g. jump in rayt is roughly 2 * Th * nbdry
    cnjump = -glm::dot(gradc, rayn_tilde - rayn);
    csjump = -glm::dot(gradc, rayt_tilde - rayt);
    
    if(isTop){
        cnjump = -cnjump; // this is because the (t,n) system of the top boundary has a different sense to the bottom boundary
        rn = -rn;
    }
    
    rm = Tg / Th; // this is tan( alpha ) where alpha is the angle of incidence
    rn = rn + rm * (RC(2.0) * cnjump - rm * csjump) / SQ(ccpx.real());
    
    if(Beam->Type[2] == 'D'){
        rn = RC(2.0) * rn;
    }else if(Beam->Type[2] == 'Z'){
        rn = RC(0.0);
    }
    
    newPoint.c   = ccpx.real();
    newPoint.tau = oldPoint.tau;
    newPoint.p   = oldPoint.p + oldPoint.q * rn;
    newPoint.q   = oldPoint.q;
    
    // account for phase change
    
    if(hs.bc == 'R'){ // rigid
        newPoint.Amp   = oldPoint.Amp;
        newPoint.Phase = oldPoint.Phase;
    }else if(hs.bc == 'V'){ // vacuum
        newPoint.Amp   = oldPoint.Amp;
        newPoint.Phase = oldPoint.Phase + M_PI;
    }else if(hs.bc == 'F'){ // file
        RInt.theta = RadDeg * STD::abs(STD::atan2(Th, Tg)); // angle of incidence (relative to normal to bathymetry)
        if(RInt.theta > RC(90.0)) RInt.theta = RC(180.0) - RInt.theta; // reflection coefficient is symmetric about 90 degrees
        InterpolateReflectionCoefficient(RInt, RefC, Npts);
        newPoint.Amp   = oldPoint.Amp * RInt.r;
        newPoint.Phase = oldPoint.Phase + RInt.phi;
    }else if(hs.bc == 'A' || hs.bc == 'G'){ // half-space
        kx = omega * Tg; // wavenumber in direction parallel      to bathymetry
        kz = omega * Th; // wavenumber in direction perpendicular to bathymetry (in ocean)
        cpx kx2 = SQ(kx);
        
        // notation below is a bit mis-leading
        // kzS, kzP is really what I called gamma in other codes, and differs by a factor of +/- i
        if(hs.cS.real() > RC(0.0)){
            kzS2 = kx2 - SQ(omega / hs.cS);
            kzP2 = kx2 - SQ(omega / hs.cP);
            kzS  = STD::sqrt(kzS2);
            kzP  = STD::sqrt(kzP2);
            mu   = hs.rho * SQ(hs.cS);
            
            y2 = (SQ(kzS2 + kx2) - RC(4.0) * kzS * kzP * kx2) * mu;
            y4 = kzP * (kx2 - kzS2);
            
            f = SQ(omega) * y4;
            g = y2;
        }else{
            kzP = STD::sqrt(kx2 - SQ(omega / hs.cP));
            
            // Intel and GFortran compilers return different branches of the SQRT for negative reals
            // LP: looks like this just means we want the positive branch
            if(kzP.real() == RC(0.0) && kzP.imag() < RC(0.0)) kzP = -kzP;
            f = kzP;
            g = hs.rho;
        }
        
        Refl = -(rho * f - J * kz * g) / (rho * f + J * kz * g); // complex reflection coef.
        
        if(STD::abs(Refl) < RC(1.0e-5)){ // kill a ray that has lost its energy in reflection
            newPoint.Amp   = RC(0.0);
            newPoint.Phase = oldPoint.Phase;
        }else{
            newPoint.Amp   = STD::abs(Refl) * oldPoint.Amp;
            newPoint.Phase = oldPoint.Phase + STD::atan2(Refl.imag(), Refl.real());
            
            // compute beam-displacement Tindle, Eq. (14)
            // needs a correction to beam-width as well ...
            // LP: most of these variables don't exist, likely very old code
            // if(kz2Sq.real() < RC(0.0)){
            //     rhoW = RC(1.0); // density of water
            //     rhoWSq = rhoW * rhoW;
            //     rhoHSSq = rhoHS * rhoHS;
            //     delta = RC(2.0) * gk * rhoW * rhoS * (kz1Sq - kz2Sq) /
            //         (kz1 * i * kz2 *
            //             (-rhoWSq * kz2Sq + rhoHSSq * kz1Sq));
            //     rv[is+1] = rv[is+1] + delta;
            // }
            
            if(Beam->Type[3] == 'S'){ // beam displacement & width change (Seongil's version)
                ch = oldPoint.c / STD::conj(hs.cP);
                co = oldPoint.t.x * oldPoint.c;
                si = oldPoint.t.y * oldPoint.c;
                ck = omega / oldPoint.c;
                
                a   = RC(2.0) * hs.rho * (RC(1.0) - SQ(ch));
                b   = SQ(co) - SQ(ch);
                d   = SQ(hs.rho) * SQ(si) + b;
                sb  = STD::sqrt(b);
                cco = SQ(co);
                ssi = SQ(si);
                
                if(si != RC(0.0)){
                    delta = a * co / si / (ck * sb * d); // Do we need an abs() on this???
                }else{
                    delta = RC(0.0);
                }
                
                pdelta = delta.real() / (oldPoint.c / co);
                // LP: The spacing in the original version of this formula,
                // the fact that several terms could be factored out to reduce
                // computation, and the repeated divisons, lead me to believe
                // that it may not be correct.
                // Here is the original version with the weird spacing:
                // ddelta = -a / (ck*sb*d) - a*cco / ssi / (ck*sb*d) + a*cco / (ck*b*sb*d)
                //     -a*co / si / (ck*sb*d*d) * (RC(2.0)* SQ(hs.rho) *si*co-RC(2.0)*co*si);
                // Here is a version with things factored better:
                cpx cksbd = ck * sb * d;
                ddelta = a * (cco / (cksbd * b)
                    - (RC(1.0) + (cco / ssi)) / cksbd
                    - RC(2.0) * SQ(co) * (SQ(hs.rho) - RC(1.0)) / (cksbd * d) );
                rddelta = -delta.real();
                sddelta = rddelta / STD::abs(rddelta);
                
                // next 3 lines have an update by Diana McCammon to allow a sloping bottom
                // I think the formulas are good, but this won't be reliable because it doesn't have the logic
                // that tracks crossing into new segments after the ray displacement.
                
                theta_bot = STD::atan(tBdry.y / tBdry.x); // bottom angle
                newPoint.x.x = newPoint.x.x + delta.real() * STD::cos(theta_bot); // range displacement
                newPoint.x.y = newPoint.x.y + delta.real() * STD::sin(theta_bot); // depth displacement
                newPoint.tau = newPoint.tau + pdelta; // phase change
                newPoint.q   = newPoint.q + sddelta * rddelta * si * ccpx.real() * oldPoint.p; // beam-width change
            }
        }
    }else{
        printf("Reflect2D: Unknown boundary condition type\n");
        bail();
    }
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
 * IsegTop, IsegBot: indices that point to the current active segment
 * rTopSeg, rBotSeg: range intervals defining the current active segment
 */
HOST_DEVICE inline bool RayInit(int32_t isrc, int32_t ialpha, real &SrcDeclAngle,
    ray2DPt &point0, vec2 &gradc, real &DistBegTop, real &DistBegBot, 
    int32_t &IsegTop, int32_t &IsegBot, vec2 &rTopSeg, vec2 &rBotSeg,
    int32_t &iSegz, int32_t &iSegr, BdryType &Bdry,
    const BdryType *ConstBdry, const BdryInfo *bdinfo, const ReflectionInfo *refl,
    const SSPStructure *ssp, const Position *Pos, const AnglesStructure *Angles,
    const FreqInfo *freqinfo, const BeamStructure *Beam, const BeamInfo *beaminfo)
{
    // LP: This part from BellhopCore
    
    float omega = RC(2.0) * M_PI * freqinfo->freq0;
    vec2 xs = vec2(RC(0.0), Pos->Sz[isrc]); // x-y coordinate of the source
    
    /*
    // LP: BUG: BELLHOP does not reinitialize these between rays, so their
    // initial state is the final state from the previous ray. Normally, the
    // initial state should not really matter, but in many examples, sources
    // are exactly on a boundary. Due to some greater-than-or-equal-to
    // versus greater-than signs, BELLHOP only searches for a new segment when
    // the point is strictly outside the current segment, but when searching
    // it uses a half-open interval. This means that if the initial state is
    // 0, it will stay with segment 0 for the first step, but if the initial
    // state is nonzero, it will move to segment 1 for the first step. This 
    // means there may or may not be one extra very small step at the beginning,
    // depending on what happened in previous rays.
    // Normally, this extra step would change very little about the results,
    // but in TL SGB mode, there is *another* bug (acknowledged by mbp) that
    // all steps are assumed to be of the default step size. Thus, adding one
    // extra very small step actually means adding one full-size step, which
    // substantially changes the results.
    // In bellhopcxx / bellhopcuda, we *cannot* have one ray's initial state be
    // dependent on another's final state, because they are processed in 
    // parallel (besides this being physically absurd). So instead, we assume
    // that ray 0 has properly initialized states, but all later rays get
    // initial states which are different from their correct segment, and so
    // they move and use the half-open search. This still does not match 100%
    // though, since hypothetically one ray could end in the same segment as the
    // next one starts, leading to the fully closed interval being used in
    // BELLHOP and possibly edge cases being different.
    if((isrc == 0 && ialpha == 0) || ssp->NPts <= 1){
        iSegz = 0;
    }else{
        // This will always be the wrong initial state (and therefore cause a
        // move using the half-open interval, like BELLHOP) unless xs.y == 
        // ssp->z[1], in which case it will be 1 but also not cause a move and
        // therefore stay 1. This should always produce results matching BELLHOP
        // as long as the previous ray really did end in some higher segment.
        iSegz = (xs.y <= ssp->z[1]) ? 1 : 0;
    }
    if((isrc == 0 && ialpha == 0) || ssp->Nr <= 1){
        iSegr = 0;
    }else{
        // Same thing except xs.x is always 0.
        iSegr = (xs.x <= ssp->Seg.r[1]) ? 1 : 0;
    }
    */
    // LP: Changed in BELLHOP to just reinitialize before every ray.
    iSegz = iSegr = 0;
    cpx ccpx; real crr, crz, czz, rho;
    EvaluateSSP(xs, ccpx, gradc, crr, crz, czz, rho, freqinfo->freq0, ssp, iSegz, iSegr);
    
    // Are there enough beams?
    real DalphaOpt = STD::sqrt(ccpx.real() / (RC(6.0) * freqinfo->freq0 * Pos->Rr[Pos->NRr-1]));
    int32_t NalphaOpt = 2 + (int)((Angles->alpha[Angles->Nalpha-1] - Angles->alpha[0]) / DalphaOpt);
    
    if(Beam->RunType[0] == 'C' && Angles->Nalpha < NalphaOpt && ialpha == 0){
        printf("Warning in bellhopcuda : Too few beams\nNalpha should be at least = %d\n", NalphaOpt);
    }
    
    real alpha = Angles->alpha[ialpha]; // initial angle
    SrcDeclAngle = RadDeg * alpha; // take-off declination angle in degrees
    
    int32_t ibp = BinarySearchLEQ(beaminfo->SrcBmPat, beaminfo->NSBPPts, 2, 0, SrcDeclAngle);
    ibp = STD::min(ibp, beaminfo->NSBPPts-2); // don't go past end of table
    
    // linear interpolation to get amplitude
    real s = (SrcDeclAngle - beaminfo->SrcBmPat[2*ibp+0]) / (beaminfo->SrcBmPat[2*(ibp+1)+0] - beaminfo->SrcBmPat[2*ibp+0]);
    float Amp0 = (RC(1.0) - s) * beaminfo->SrcBmPat[2*ibp+1] + s * beaminfo->SrcBmPat[2*(ibp+1)+1]; // initial amplitude
    
    // Lloyd mirror pattern for semi-coherent option
    if(Beam->RunType[0] == 'S')
        Amp0 *= STD::sqrt(RC(2.0)) * STD::abs(STD::sin(omega / ccpx.real() * xs.y * STD::sin(alpha)));
        
    // LP: This part from TraceRay2D
    
    point0 = ray2DPt{
        /*.NumTopBnc =*/ 0,
        /*.NumBotBnc =*/ 0,
        /*.x         =*/ xs,
        /*.t         =*/ vec2(STD::cos(alpha), STD::sin(alpha)) / ccpx.real(),
        /*.p         =*/ vec2(RC(1.0), RC(0.0)),
        /*.q         =*/ vec2(RC(0.0), RC(1.0)),
        /*.c         =*/ ccpx.real(),
        /*.Amp       =*/ Amp0,
        /*.Phase     =*/ RC(0.0),
        /*.tau       =*/ cpx(RC(0.0), RC(0.0)),
    };
    // second component of qv is not used in geometric beam tracing
    // set I.C. to 0 in hopes of saving run time
    if(Beam->RunType[1] == 'G') point0.q = vec2(RC(0.0), RC(0.0));
    
    IsegTop = 0; IsegBot = 0;
    GetTopSeg(xs.x, IsegTop, rTopSeg, bdinfo); // identify the top    segment above the source
    GetBotSeg(xs.x, IsegBot, rBotSeg, bdinfo); // identify the bottom segment below the source
    
    // convert range-dependent geoacoustic parameters from user to program units
    // LP: BELLHOP uses all values from ConstBdry except replaces cP, cS, and rho
    // from the current segment in bdinfo. rho is read from files in ReadATI and
    // Read BTY, and cP and cS are computed in core_setup. bc, which is also read
    // by Reflect2D, is never set in bdinfo and is left alone from ConstBdry.
    Bdry = *ConstBdry;
    if(bdinfo->atiType[1] == 'L') CopyHSInfo(Bdry.Top.hs, bdinfo->Top[IsegTop].hs);
    if(bdinfo->btyType[1] == 'L') CopyHSInfo(Bdry.Bot.hs, bdinfo->Bot[IsegBot].hs);
    
    vec2 dEndTop_temp, dEndBot_temp;
    Distances2D(point0.x, bdinfo->Top[IsegTop].x, bdinfo->Bot[IsegBot].x, 
        dEndTop_temp, dEndBot_temp,
        bdinfo->Top[IsegTop].n, bdinfo->Bot[IsegBot].n, DistBegTop, DistBegBot);
        
    if(DistBegTop <= RC(0.0) || DistBegBot <= RC(0.0)){
        printf("Terminating the ray trace because the source is on or outside the boundaries\n");
        // printf("xs (%g,%g) Bot.x (%g,%g) Bot.n (%g,%g) DistBegBot %g\n",
        //     xs.x, xs.y,
        //     bdinfo->Bot[IsegBot].x.x, bdinfo->Bot[IsegBot].x.y,
        //     bdinfo->Bot[IsegBot].n.x, bdinfo->Bot[IsegBot].n.y, DistBegBot);
        return false; // source must be within the medium
    }
    
    return true;
}

/**
 * LP: Pulled out contents of ray update loop. Returns the number of ray steps
 * taken (i.e. normally 1, or 2 if reflected).
 */
HOST_DEVICE inline int32_t RayUpdate(
    const ray2DPt &point0, ray2DPt &point1, ray2DPt &point2,
    const real &DistBegTop, const real &DistBegBot, real &DistEndTop, real &DistEndBot, 
    int32_t &IsegTop, int32_t &IsegBot, vec2 &rTopSeg, vec2 &rBotSeg,
    int32_t &iSmallStepCtr, int32_t &iSegz, int32_t &iSegr,
    BdryType &Bdry, const BdryInfo *bdinfo, const ReflectionInfo *refl,
    const SSPStructure *ssp, const FreqInfo *freqinfo, const BeamStructure *Beam)
{
    int32_t numRaySteps = 1;
    
    Step2D(point0, point1, bdinfo->Top[IsegTop].x, bdinfo->Top[IsegTop].n,
        bdinfo->Bot[IsegBot].x, bdinfo->Bot[IsegBot].n, rTopSeg, rBotSeg, 
        freqinfo->freq0, Beam, ssp, iSegz, iSegr, iSmallStepCtr);
    
    // New altimetry segment?
    if(point1.x.x < rTopSeg.x || point1.x.x > rTopSeg.y){
        GetTopSeg(point1.x.x, IsegTop, rTopSeg, bdinfo);
        if(bdinfo->atiType[1] == 'L') CopyHSInfo(Bdry.Top.hs, bdinfo->Top[IsegTop].hs); // grab the geoacoustic info for the new segment
    }
    
    // New bathymetry segment?
    if(point1.x.x < rBotSeg.x || point1.x.x > rBotSeg.y){
        GetBotSeg(point1.x.x, IsegBot, rBotSeg, bdinfo);
        if(bdinfo->btyType[1] == 'L') CopyHSInfo(Bdry.Bot.hs, bdinfo->Bot[IsegBot].hs); // grab the geoacoustic info for the new segment
    }
    
    // Reflections?
    // Tests that ray at step is is inside, and ray at step is+1 is outside
    // to detect only a crossing from inside to outside
    // DistBeg is the distance at point0, which is saved
    // DistEnd is the distance at point1, which needs to be calculated
    vec2 dEndTop, dEndBot;
    Distances2D(point1.x, bdinfo->Top[IsegTop].x, bdinfo->Bot[IsegBot].x,
        dEndTop, dEndBot,
        bdinfo->Top[IsegTop].n, bdinfo->Bot[IsegBot].n, DistEndTop, DistEndBot);
    
    bool topRefl = DistBegTop > RC(0.0) && DistEndTop <= RC(0.0); // test top reflection
    bool botRefl = DistBegBot > RC(0.0) && DistEndBot <= RC(0.0); // test bottom reflection
    // LP: Merging these cases is important for GPU performance.
    if(topRefl || botRefl){
        vec2 dEnd = topRefl ? dEndTop : dEndBot;
        BdryPtFull *bd0 = topRefl ? &bdinfo->Top[IsegTop] : &bdinfo->Bot[IsegBot];
        BdryPtFull *bd1 = &bd0[1]; // LP: next segment
        vec2 nInt, tInt;
        // LP: FORTRAN actually checks if the whole string is just "C", not just the first char
        if((topRefl ? bdinfo->atiType[0] : bdinfo->btyType[0]) == 'C'){
            real sss = glm::dot(dEnd, bd0->t) / bd0->Len; // proportional distance along segment
            nInt = (RC(1.0) - sss) * bd0->Noden + sss * bd1->Noden;
            tInt = (RC(1.0) - sss) * bd0->Nodet + sss * bd1->Nodet;
        }else{
            nInt = bd0->n; // normal is constant in a segment
            tInt = bd0->t;
        }
        Reflect2D(point1, point2, 
            topRefl ? Bdry.Top.hs : Bdry.Bot.hs,
            topRefl, tInt, nInt, bd0->kappa, freqinfo->freq0,
            topRefl ? refl->RTop : refl->RBot,
            topRefl ? refl->NTopPts : refl->NBotPts,
            Beam, ssp, iSegz, iSegr);
        //Incrementing bounce count moved to Reflect2D
        numRaySteps = 2;
        Distances2D(point2.x, bdinfo->Top[IsegTop].x, bdinfo->Bot[IsegBot].x, dEndTop, dEndBot,
            bdinfo->Top[IsegTop].n, bdinfo->Bot[IsegBot].n, DistEndTop, DistEndBot);
    }
    
    return numRaySteps;
}

/**
 * Has the ray left the box, lost its energy, escaped the boundaries, or 
 * exceeded storage limit?
 * LP: Also updates DistBegTop, DistBegBot.
 */
HOST_DEVICE inline bool RayTerminate(const ray2DPt &point, int32_t &Nsteps, int32_t is,
    real &DistBegTop, real &DistBegBot, const real &DistEndTop, const real &DistEndBot,
    const BeamStructure *Beam
    )
{
    bool leftbox = STD::abs(point.x.x) > Beam->Box.r ||
                   STD::abs(point.x.y) > Beam->Box.z;
    bool lostenergy = point.Amp < RC(0.005);
    bool escapedboundaries = (DistBegTop < RC(0.0) && DistEndTop < RC(0.0)) ||
                             (DistBegBot < RC(0.0) && DistEndBot < RC(0.0));
    //bool backward = ray2D[is+1].t.x < 0; // this last test kills off a backward traveling ray
    if(leftbox || lostenergy || escapedboundaries // || backward
    ){
        /*
        if(leftbox){
            printf("Ray left beam box (%g,%g)\n", Beam->Box.r, Beam->Box.z);
        }else if(lostenergy){
            printf("Ray energy dropped to %g\n", ray2D[is+1].Amp);
        }else{
            printf("Ray escaped boundaries DistBegTop %g DistEndTop %g DistBegBot %g DistEndBot %g\n",
                DistBegTop, DistEndTop, DistBegBot, DistEndBot);
        }
        */
        Nsteps = is + 1;
        return true;
    }else if(is >= MaxN - 3){
        printf("Warning in TraceRay2D: Insufficient storage for ray trajectory\n");
        Nsteps = is;
        return true;
    }
    
    DistBegTop = DistEndTop;
    DistBegBot = DistEndBot;
    return false;
}
