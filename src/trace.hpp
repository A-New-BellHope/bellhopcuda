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
HOST_DEVICE inline void Reflect2D(int32_t &is, const HSInfo &hs, bool isTop,
    const vec2 &tBdry, const vec2 &nBdry, real kappa, real freq,
    const ReflectionCoef *RefC, int32_t Npts,
    ray2DPt *ray2D, const BeamStructure *Beam, const SSPStructure *ssp,
    int32_t &iSegz, int32_t &iSegr)
{
    int32_t is1;
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
    
    is  = is + 1;
    is1 = is + 1;
    
    Tg = glm::dot(ray2D[is].t, tBdry); // component of ray tangent, along boundary
    Th = glm::dot(ray2D[is].t, nBdry); // component of ray tangent, normal to boundary
    
    ray2D[is1].NumTopBnc = ray2D[is].NumTopBnc;
    ray2D[is1].NumBotBnc = ray2D[is].NumBotBnc;
    ray2D[is1].x         = ray2D[is].x;
    ray2D[is1].t         = ray2D[is].t - RC(2.0) * Th * nBdry; // changing the ray direction
    
    // Calculate the change in curvature
    // Based on formulas given by Muller, Geoph. J. R.A.S., 79 (1984).
    
    EvaluateSSP(ray2D[is].x, ccpx, gradc, crr, crz, czz, rho, freq, ssp, iSegz, iSegr); // just to get c [LP: ccpx.real()]
    
    // incident unit ray tangent and normal
    rayt = ccpx.real() * ray2D[is].t; // unit tangent to ray
    rayn = vec2(-rayt.y, rayt.x);     // unit normal  to ray
    
    // reflected unit ray tangent and normal (the reflected tangent, normal system has a different orientation)
    rayt_tilde = ccpx.real() * ray2D[is1].t;         // unit tangent to ray
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
    
    ray2D[is1].c   = ccpx.real();
    ray2D[is1].tau = ray2D[is].tau;
    ray2D[is1].p   = ray2D[is].p + ray2D[is].q * rn;
    ray2D[is1].q   = ray2D[is].q;
    
    // account for phase change
    
    if(hs.bc == 'R'){ // rigid
        ray2D[is1].Amp   = ray2D[is].Amp;
        ray2D[is1].Phase = ray2D[is].Phase;
    }else if(hs.bc == 'V'){ // vacuum
        ray2D[is1].Amp   = ray2D[is].Amp;
        ray2D[is1].Phase = ray2D[is].Phase + M_PI;
    }else if(hs.bc == 'F'){ // file
        RInt.theta = RadDeg * STD::abs(STD::atan2(Th, Tg)); // angle of incidence (relative to normal to bathymetry)
        if(RInt.theta > RC(90.0)) RInt.theta = RC(180.0) - RInt.theta; // reflection coefficient is symmetric about 90 degrees
        InterpolateReflectionCoefficient(RInt, RefC, Npts);
        ray2D[is1].Amp   = ray2D[is].Amp * RInt.r;
        ray2D[is1].Phase = ray2D[is].Phase + RInt.phi;
    }else if(hs.bc == 'A' || hs.bc == 'G'){ // half-space
        kx = omega * Tg; // wavenumber in direction parallel      to bathymetry
        kz = omega * Th; // wavenumber in direction perpendicular to bathymetry (in ocean)
        real kx2 = SQ(kx);
        
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
        
        Refl = -cpx(rho * f, -(kz * g)) / cpx(rho * f, kz * g); // complex reflection coef.
        
        if(STD::abs(Refl) < RC(1.0e-5)){ // kill a ray that has lost its energy in reflection
            ray2D[is1].Amp   = RC(0.0);
            ray2D[is1].Phase = ray2D[is].Phase;
        }else{
            ray2D[is1].Amp   = STD::abs(Refl) * ray2D[is].Amp;
            ray2D[is1].Phase = ray2D[is].Phase + STD::atan2(Refl.imag(), Refl.real());
            
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
                ch = ray2D[is].c / STD::conj(hs.cP);
                co = ray2D[is].t.x * ray2D[is].c;
                si = ray2D[is].t.y * ray2D[is].c;
                ck = omega / ray2D[is].c;
                
                a   = RC(2.0) * hs.rho * (RC(1.0) - SQ(ch));
                b   = SQ(co) - SQ(ch);
                d   = SQ(hs.rho) * SQ(si) + b;
                sb  = STD::sqrt(b);
                cco = SQ(co);
                ssi = SQ(si);
                
                if(si != RC(0.0)){
                    delta = a * co / si / (ck * sb * d); // Do we need an abs() on this???
                }else{
                    delta = 0.0
                }
                
                pdelta = delta.real() / (ray2D[is].c / co);
                // LP: The spacing in the original version of this formula,
                // the fact that several terms could be factored out to reduce
                // computation, and the repeated divisons, lead me to believe
                // that it may not be correct.
                // Here is the original version with the weird spacing:
                // ddelta = -a / (ck*sb*d) - a*cco / ssi / (ck*sb*d) + a*cco / (ck*b*sb*d)
                //     -a*co / si / (ck*sb*d*d) * (RC(2.0)* SQ(hs.rho) *si*co-RC(2.0)*co*si);
                // Here is a version with things factored better:
                real cksbd = ck * sb * d;
                ddelta = a * (cco / (cksbd * b)
                    - (RC(1.0) + (cco / ssi)) / cksbd
                    - RC(2.0) * SQ(co) * (SQ(hs.rho) - RC(1.0)) / (cksbd * d) );
                rddelta = -delta.real();
                sddelta = rddelta / STD::abs(rddelta);
                
                // next 3 lines have an update by Diana McCammon to allow a sloping bottom
                // I think the formulas are good, but this won't be reliable because it doesn't have the logic
                // that tracks crossing into new segments after the ray displacement.
                
                theta_bot = STD::atan(tBdry.y / tBdry.x); // bottom angle
                ray2D[is1].x.x = ray2D[is1].x.x + delta.real() * STD::cos(theta_bot); // range displacement
                ray2D[is1].x.y = ray2D[is1].x.y + delta.real() * STD::sin(theta_bot); // depth displacement
                ray2D[is1].tau = ray2D[is1].tau + pdelta; // phase change
                ray2D[is1].q   = ray2D[is1].q + sddelta * rddelta * si * c * ray2D[is].p; // beam-width change
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
 * Traces the beam corresponding to a particular take-off angle
 * 
 * xs: x-y coordinate of the source
 * alpha: initial angle
 * Amp0: initial amplitude
 */
HOST_DEVICE inline void TraceRay2D(vec2 xs, real alpha, real Amp0,
    ray2DPt *ray2D, int32_t &Nsteps,
    const BdryType *ConstBdry, const BdryInfo *bdinfo, const ReflectionInfo *refl,
    const SSPStructure *ssp, const FreqInfo *freqinfo, const BeamStructure *Beam)
{
    int32_t is, is1; // index for a step along the ray
    cpx ccpx;
    vec2 gradc, dEndTop, dEndBot, TopnInt, BotnInt, ToptInt, BottInt;
    real crr, crz, czz, rho;
    real DistBegTop, DistEndTop, DistBegBot, DistEndBot; // Distances from ray beginning, end to top and bottom
    real sss;
    
    int32_t IsegTop, IsegBot; // indices that point to the current active segment
    vec2 rTopseg, rBotseg; // range intervals defining the current active segment
    
    // Initial conditions
    
    int32_t iSmallStepCtr = 0, iSegz = 0, iSegr = 0;
    EvaluateSSP(xs, ccpx, gradc, crr, crz, czz, rho, freqinfo->freq0, iSegz, iSegr);
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
    // LP: BELLHOP uses all values from ConstBdry except replaces cP, cS, and rho
    // from the current segment in bdinfo. rho is read from files in ReadATI and
    // Read BTY, and cP and cS are computed in core_setup. bc, which is also read
    // by Reflect2D, is never set in bdinfo and is left alone from ConstBdry.
    BdryType Bdry = *ConstBdry;
    if(bdinfo->atiType[1] == 'L') CopyHSInfo(Bdry.Top.hs, bdinfo->Top[IsegTop].hs);
    if(bdinfo->btyType[1] == 'L') CopyHSInfo(Bdry.Bot.hs, bdinfo->Bot[IsegBot].hs);
    
    // Trace the beam (note that Reflect alters the step index is)
    is = 0;
    Distances2D(ray2D[0].x, bdinfo->Top[IsegTop].x, bdinfo->Bot[IsegBot].x, dEndTop, dEndBot,
        bdinfo->Top[IsegTop].n, bdinfo->Bot[IsegBot].n, DistBegTop, DistBegBot);
    
    if(DistBegTop <= 0 || DistBegBot <= 0){
        Nsteps = 1;
        printf("Terminating the ray trace because the source is on or outside the boundaries\n");
        return; // source must be within the medium
    }
    
    for(int32_t istep = 0; istep<MaxN-1; ++istep){
        is  = istep + 1;
        is1 = istep + 1;
        
        Step2D(ray2D[is], &ray2D[is1], bdinfo->Top[IsegTop].x, bdinfo->Top[IsegTop].n,
            bdinfo->Bot[IsegBot].x, bdinfo->Bot[IsegBot].n, freqinfo->freq0, Beam, ssp, 
            iSegz, iSegr, iSmallStepCtr);
        
        // New altimetry segment?
        if(ray2D[is1].x.x < rTopSeg.x || ray2D[is1].x.x > rTopSeg.y){
            GetTopSeg(ray2D[is1].x.x, IsegTop, rTopseg);
            if(bdinfo->atiType[1] == 'L') CopyHSInfo(Bdry.Top.hs, bdinfo->Top[IsegTop].hs); // grab the geoacoustic info for the new segment
        }
        
        // New bathymetry segment?
        if(ray2D[is1].x.x < rBotSeg.x || ray2D[is1].x.x > rBotSeg.y){
            GetBotSeg(ray2D[is1].x.x, IsegBot, rBotseg);
            if(bdinfo->btyType[1] == 'L') CopyHSInfo(Bdry.Bot.hs, bdinfo->Bot[IsegBot].hs); // grab the geoacoustic info for the new segment
        }
        
        // Reflections?
        // Tests that ray at step is is inside, and ray at step is+1 is outside
        // to detect only a crossing from inside to outside
        // DistBeg is the distance at step is,   which is saved
        // DistEnd is the distance at step is+1, which needs to be calculated
        
        Distances2D(ray2D[is1].x, bdinfo->Top[IsegTop].x, bdinfo->Bot[IsegBot].x,
            dEndTop, dEndBot,
            bdinfo->Top[IsegTop].n, bdinfo->Bot[IsegBot].n, DistEndTop, DistEndBot);
        
        if(DistBegTop > RC(0.0) && DistEndTop <= RC(0.0)){ // test top reflection
            
            if(bdinfo->atiType[0] == 'C'){ // LP: Actually checking if the whole string is just "C", not just the first char
                sss = glm::dot(dEndTop, bdinfo->Top[IsegTop].t) / bdinfo->Top[IsegTop].Len; // proportional distance along segment
                TopnInt = (RC(1.0) - sss) * bdinfo->Top[IsegTop].Noden + sss * bdinfo->Top[IsegTop+1].Noden;
                ToptInt = (RC(1.0) - sss) * bdinfo->Top[IsegTop].Nodet + sss * bdinfo->Top[IsegTop+1].Nodet;
            }else{
                TopnInt = bdinfo->Top[IsegTop].n; // normal is constant in a segment
                ToptInt = bdinfo->Top[IsegTop].t;
            }
            
            Reflect2D(is, Bdry.Top.hs, true, ToptInt, TopnInt, bdinfo->Top[IsegTop].Kappa, 
                freqinfo->freq0, refl->RTop, refl->NTopPTS, ray2D, Beam, ssp, iSegz, iSegr);
            ++ray2D[is+1].NumTopBnc;
            
            Distances2D(ray2D[is+1].x, bdinfo->Top[IsegTop].x, bdinfo->Bot[IsegBot].x, dEndTop, dEndBot,
                bdinfo->Top[IsegTop].n, bdinfo->Bot[IsegBot].n, DistEndTop, DistEndBot);
                
        }else if(DistBegBot > RC(0.0) && DistEndBot <= RC(0.0)){ // test bottom reflection
            
            if(bdinfo->btyType[0] == 'C'){ // LP: Actually checking if the whole string is just "C", not just the first char
                sss = glm::dot(dEndBot, bdinfo->Bot[IsegBot].t) / bdinfo->Bot[IsegBot].Len; // proportional distance along segment
                BotnInt = (RC(1.0) - sss) * bdinfo->Bot[IsegBot].Noden + sss * bdinfo->Bot[IsegBot+1].Noden;
                BottInt = (RC(1.0) - sss) * bdinfo->Bot[IsegBot].Nodet + sss * bdinfo->Bot[IsegBot+1].Nodet;
            }else{
                BotnInt = bdinfo->Bot[IsegBot].n; // normal is constant in a segment
                BottInt = bdinfo->Bot[IsegBot].t;
            }
            
            Reflect2D(is, Bdry.Bot.hs, false, BottInt, BotnInt, bdinfo->Bot[IsegBot].Kappa, 
                freqinfo->freq0, refl->RBot, refl->NBotPTS, ray2D, Beam, ssp, iSegz, iSegr);
            ++ray2D[is+1].NumBotBnc;
            
            Distances2D(ray2D[is+1].x, bdinfo->Top[IsegTop].x, bdinfo->Bot[IsegBot].x, dEndTop, dEndBot,
                bdinfo->Top[IsegTop].n, bdinfo->Bot[IsegBot].n, DistEndTop, DistEndBot);
            
        }
        
        // Has the ray left the box, lost its energy, escaped the boundaries, or exceeded storage limit?
        if( STD::abs(ray2D[is+1].x.x) > Beam->Box.r ||
            STD::abs(ray2d[is+1].x.y) > Beam->Box.z ||
            ray2d[is+1].Amp < RC(0.005) ||
            (DistBegTop < RC(0.0) && DistEndTop < RC(0.0)) ||
            (DistBegBot < RC(0.0) && DistEndBot < RC(0.0)) 
            // ray2d[is+1].t.x < 0 // this last test kills off a backward traveling ray
        ){
            Nsteps = is + 1;
            break;
        }else if(is >= MaxN - 4){
            printf("Warning in TraceRay2D: Insufficient storage for ray trajectory\n");
            Nsteps = is;
            break;
        }
        
        DistBegTop = DistEndTop;
        DistBegBot = DistEndBot;        
    }
}

//TODO 
/*
switch(Beam->RunType[0]){
case 'C':
case 'S':
case 'I':
    // TL calculation, zero out pressure matrix
    memset(inflinfo->u, 0, TODO);
    break;
case 'A':
case 'a':
    // Arrivals calculation, zero out arrival matrix
    memset(arrinfo->Narr, 0, TODO);
    break;
}
*/



//TODO only call for valid ialpha w.r.t. Angles->iSingleAlpha
inline void CoreSingleBeam(int32_t is, int32_t ialpha, ray2DPt *ray2D,
    real &SrcDeclAngle, int32_t &Nsteps,
    const BdryType *Bdry, const BdryInfo *bdinfo, const ReflectionInfo *refl,
    const SSPStructure *ssp, const Position *Pos, const AnglesStructure *Angles,
    const FreqInfo *freqinfo, const BeamStructure *Beam, const BeamInfo *beaminfo)
{
    float omega = RC(2.0) * M_PI * freqinfo->freq0;
    
    vec2 xs = vec2(RC(0.0), Pos->sz[is]);
    
    cpx ccpx;
    EvaluateSSP(xs, ccpx, freqinfo->freq0, ssp, iSegz, iSegr);
    real RadMax = RC(5.0) * ccpx.real() / freqinfo->freq0; // 5 wavelength max radius
    
    // Are there enough beams?
    real DalphaOpt = STD::sqrt(ccpx.real() / (RC(6.0) * freqinfo->freq0 * Pos->Rr[Pos-NRr-1]);
    int32_t NalphaOpt = 2 + (int)((Angles->alpha[Angles->Nalpha-1] - Angles->alpha[0]) / DalphaOpt);
    
    if(Beam->RunType[0] == 'C' && Angles->Nalpha < NalphaOpt){
        printf("Warning in bellhopcuda : Too few beams\nNalpha should be at least = %d\n", NalphaOpt);
    }
    
    real alpha = Angles->alpha[ialpha];
    SrcDeclAngle = RadDeg * alpha; // take-off declination angle in degrees
    
    int32_t ibp = BinarySearch(beaminfo->SrcBmPat, beaminfo->NSBPPts, 2, 0, SrcDeclAngle);
    ibp = STD::min(ibp, beaminfo->NSBPPts-2); // don't go past end of table
    
    // linear interpolation to get amplitude
    real s = (SrcDeclAngle - beaminfo->SrcBmPat[2*ibp+0]) / (beaminfo->SrcBmPat[2*(ibp+1)+0] - beaminfo->SrcBmPat[2*ibp+0]);
    float Amp0 = (RC(1.0) - s) * beaminfo->SrcBmPat[2*ibp+1] + s * beaminfo->SrcBmPat[2*(ibp+1)+1];
    
    // Lloyd mirror pattern for semi-coherent option
    if(Beam->RunType[0] == 'S')
        Amp0 *= STD::sqrt(RC(2.0)) * STD::abs(STD::sin(omega / ccpx.real() * xs.y * STD::sin(alpha)));
        
    // show progress ...
    //LP: Removed
    // if(((ialpha - 1) % STD::max(Angles->Nalpha / 50, 1)) == 0){
    //     printf("Tracing beam %d %f\n", ialpha, SrcDeclAngle);
    // }
    
    TraceRay2D(xs, alpha, Amp0, ray2D, Nsteps, Bdry, bdinfo, refl, ssp, freqinfo, Beam);
}
