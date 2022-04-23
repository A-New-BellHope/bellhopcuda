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

namespace bhc {

/**
 * Given an angle RInt%ThetaInt, returns the magnitude and
 * phase of the reflection coefficient (RInt%R, RInt%phi).
 * 
 * Uses linear interpolation using the two nearest abscissas
 * Assumes phi has been unwrapped so that it varies smoothly.
 * I tried modifying it to allow a complex angle of incidence but
 * stopped when I realized I needed to fuss with a complex atan2 routine
 * LP: C++ supports complex atan, though there is no complex atan2.
 * 
 * RInt: interpolated value of refl. coef.
 * r: Reflection coefficient table
 * NPts: # pts in refl. coef.
 */
HOST_DEVICE inline void InterpolateReflectionCoefficient(ReflectionCoef &RInt,
    const ReflectionCoef *r, int32_t NPts)
{
    int32_t iLeft, iRight, iMid;
    real alpha, thetaIntr;
    
    iLeft = 0;
    iRight = NPts - 1;
    
    // LP: This was originally the FORTRAN version of RInt.theta.real(), but
    // theta is definitely already a real (originally double).
    thetaIntr = RInt.theta; // This should be unnecessary? probably used when I was doing complex angles
    
    // Three cases: ThetaInt left, in, or right of tabulated interval
    
    if(thetaIntr < r[iLeft].theta){
        // iRight = 1;
        RInt.r   = FL(0.0); // r[iLeft].r
        RInt.phi = FL(0.0); // r[iLeft].phi
        printf("Warning in InterpolateReflectionCoefficient : Refl. Coef. being "
            "set to 0 outside tabulated domain : angle = %f, lower limit = %f",
            thetaIntr, r[iLeft].theta);
    }else if(thetaIntr > r[iRight].theta){
        // iLeft = NPts - 2;
        RInt.r   = FL(0.0); // r[iRight].r
        RInt.phi = FL(0.0); // r[iRight].phi
        // printf("Warning in InterpolateReflectionCoefficient : Refl. Coef. being "
        //     "set to 0 outside tabulated domain : angle = %f, lower limit = %f",
        //     thetaIntr, r[iRight].theta);
    }else{
        // Search for bracketing abscissas: Log2( NPts ) stabs required for a bracket
        
        while(iLeft != iRight - 1){
            iMid = (iLeft + iRight) / 2;
            if(r[iMid].theta > thetaIntr){
                iRight = iMid;
            }else{
                iLeft = iMid;
            }
        }
        
        // Linear interpolation for reflection coef
        
        alpha    = (RInt.theta - r[iLeft].theta) / (r[iRight].theta - r[iLeft].theta);
        RInt.r   = (FL(1.0) - alpha) * r[iLeft].r   + alpha * r[iRight].r;
        RInt.phi = (FL(1.0) - alpha) * r[iLeft].phi + alpha * r[iRight].phi;
    }
}

/**
 * Optionally read in reflection coefficient for Top or Bottom boundary
 * 
 * BotRC, TopRC: flag set to 'F' if refl. coef. is to be read from a File
 */
inline void ReadReflectionCoefficient(std::string FileRoot, char BotRC, char TopRC,
    PrintFileEmu &PRTFile, ReflectionInfo *refl)
{
    if(BotRC == 'F'){
        PRTFile << "__________________________________________________________________________\n\n";
        PRTFile << "Using tabulated bottom reflection coef.\n";
        LDIFile BRCFile(FileRoot + ".brc");
        if(!BRCFile.Good()){
            PRTFile << "BRCFile = " << FileRoot + ".brc\n";
            std::cout << "ReadReflectionCoefficient: Unable to open Bottom Reflection Coefficient file\n";
            std::abort();
        }
        
        LIST(BRCFile); BRCFile.Read(refl->NBotPts);
        PRTFile << "Number of points in bottom reflection coefficient = " << refl->NBotPts << "\n";
        
        checkallocate(refl->RBot, refl->NBotPts);
        
        LIST(BRCFile);
        for(int32_t itheta=0; itheta<refl->NBotPts; ++itheta){
            BRCFile.Read(refl->RBot[itheta].theta);
            BRCFile.Read(refl->RBot[itheta].r);
            BRCFile.Read(refl->RBot[itheta].phi);
            refl->RBot->phi *= DegRad; // convert to radians
        }
    }else{ // should allocate something anyway, since variable is passed
        checkallocate(refl->RBot, 1);
    }
    
    // Optionally read in top reflection coefficient
    
    if(TopRC == 'F'){
        PRTFile << "__________________________________________________________________________\n\n";
        PRTFile << "Using tabulated top    reflection coef.\n";
        LDIFile TRCFile(FileRoot + ".trc");
        if(!TRCFile.Good()){
            PRTFile << "TRCFile = " << FileRoot + ".trc\n";
            std::cout << "ReadReflectionCoefficient: Unable to open Top Reflection Coefficient file\n";
            std::abort();
        }
        
        LIST(TRCFile); TRCFile.Read(refl->NTopPts);
        PRTFile << "Number of points in top reflection coefficient = " << refl->NTopPts << "\n";
        
        checkallocate(refl->RTop, refl->NTopPts);
        
        LIST(TRCFile);
        for(int32_t itheta=0; itheta<refl->NTopPts; ++itheta){
            TRCFile.Read(refl->RTop[itheta].theta);
            TRCFile.Read(refl->RTop[itheta].r);
            TRCFile.Read(refl->RTop[itheta].phi);
            refl->RTop->phi *= DegRad; // convert to radians
        }
    }else{ // should allocate something anyway, since variable is passed
        checkallocate(refl->RTop, 1);
    }
    
    // Optionally read in internal reflection coefficient data
    
    if(BotRC == 'P'){
        std::cout << "Internal reflections not supported by BELLHOP and therefore "
            "not supported by " BHC_PROGRAMNAME "\n";
        std::abort();
    }
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
HOST_DEVICE inline void Reflect2D(const rayPt<false> &oldPoint, rayPt<false> &newPoint,
    const HSInfo &hs, bool isTop,
    const vec2 &tBdry, const vec2 &nBdry, real kappa, real freq,
    const ReflectionCoef *RefC, int32_t Npts,
    const BeamStructure *Beam, const SSPStructure *ssp, SSPSegState &iSeg)
{
    SSPOutputs<false> o;
    real rm, rn, Tg, Th;
    vec2 rayt, rayn, rayt_tilde, rayn_tilde;
    real cnjump, csjump; // for curvature change
    real ck, co, si, cco, ssi, pdelta, rddelta, sddelta, theta_bot; // for beam shift
    cpx kx, kz, kzP, kzS, kzP2, kzS2, mu, f, g, y2, y4, Refl; // for tabulated reflection coef.
    cpx ch, a, b, d, sb, delta, ddelta; // for beam shift
    ReflectionCoef RInt;
    real omega = FL(2.0) * REAL_PI * freq;
    
    Tg = glm::dot(oldPoint.t, tBdry); // component of ray tangent, along boundary
    Th = glm::dot(oldPoint.t, nBdry); // component of ray tangent, normal to boundary
    
    // LP: Incrementing bounce count moved here
    newPoint.NumTopBnc = oldPoint.NumTopBnc + (isTop ? 1 : 0);
    newPoint.NumBotBnc = oldPoint.NumBotBnc + (isTop ? 0 : 1);
    newPoint.x         = oldPoint.x;
    newPoint.t         = oldPoint.t - FL(2.0) * Th * nBdry; // changing the ray direction
    
    // Calculate the change in curvature
    // Based on formulas given by Muller, Geoph. J. R.A.S., 79 (1984).
    
    // just to get c [LP: o.ccpx.real(); also, this is wrong, it is also using o.gradc]
    EvaluateSSP<false, false>(oldPoint.x, oldPoint.t, o, ssp, iSeg);
    
    // incident unit ray tangent and normal
    rayt = o.ccpx.real() * oldPoint.t; // unit tangent to ray
    rayn = vec2(-rayt.y, rayt.x);     // unit normal  to ray
    
    // reflected unit ray tangent and normal (the reflected tangent, normal system has a different orientation)
    rayt_tilde = o.ccpx.real() * newPoint.t;         // unit tangent to ray
    rayn_tilde = -vec2(-rayt_tilde.y, rayt_tilde.x); // unit normal  to ray
    
    rn = FL(2.0) * kappa / SQ(o.ccpx.real()) / Th; // boundary curvature correction
    
    // get the jumps (this could be simplified, e.g. jump in rayt is roughly 2 * Th * nbdry
    cnjump = -glm::dot(o.gradc, rayn_tilde - rayn);
    csjump = -glm::dot(o.gradc, rayt_tilde - rayt);
    
    if(isTop){
        cnjump = -cnjump; // this is because the (t,n) system of the top boundary has a different sense to the bottom boundary
        rn = -rn;
    }
    
    rm = Tg / Th; // this is tan( alpha ) where alpha is the angle of incidence
    rn = rn + rm * (FL(2.0) * cnjump - rm * csjump) / SQ(o.ccpx.real());
    
    if(Beam->Type[2] == 'D'){
        rn = FL(2.0) * rn;
    }else if(Beam->Type[2] == 'Z'){
        rn = FL(0.0);
    }
    
    newPoint.c   = o.ccpx.real();
    newPoint.tau = oldPoint.tau;
    newPoint.p   = oldPoint.p + oldPoint.q * rn;
    newPoint.q   = oldPoint.q;
    
    // account for phase change
    
    if(hs.bc == 'R'){ // rigid
        newPoint.Amp   = oldPoint.Amp;
        newPoint.Phase = oldPoint.Phase;
    }else if(hs.bc == 'V'){ // vacuum
        newPoint.Amp   = oldPoint.Amp;
        newPoint.Phase = oldPoint.Phase + REAL_PI;
    }else if(hs.bc == 'F'){ // file
        RInt.theta = RadDeg * STD::abs(STD::atan2(Th, Tg)); // angle of incidence (relative to normal to bathymetry)
        if(RInt.theta > FL(90.0)) RInt.theta = FL(180.0) - RInt.theta; // reflection coefficient is symmetric about 90 degrees
        InterpolateReflectionCoefficient(RInt, RefC, Npts);
        newPoint.Amp   = oldPoint.Amp * RInt.r;
        newPoint.Phase = oldPoint.Phase + RInt.phi;
    }else if(hs.bc == 'A' || hs.bc == 'G'){ // half-space
        kx = omega * Tg; // wavenumber in direction parallel      to bathymetry
        kz = omega * Th; // wavenumber in direction perpendicular to bathymetry (in ocean)
        cpx kx2 = SQ(kx);
        
        // notation below is a bit mis-leading
        // kzS, kzP is really what I called gamma in other codes, and differs by a factor of +/- i
        if(hs.cS.real() > FL(0.0)){
            kzS2 = kx2 - SQ(omega / hs.cS);
            kzP2 = kx2 - SQ(omega / hs.cP);
            kzS  = STD::sqrt(kzS2);
            kzP  = STD::sqrt(kzP2);
            mu   = hs.rho * SQ(hs.cS);
            
            y2 = (SQ(kzS2 + kx2) - FL(4.0) * kzS * kzP * kx2) * mu;
            y4 = kzP * (kx2 - kzS2);
            
            f = SQ(omega) * y4;
            g = y2;
        }else{
            kzP = STD::sqrt(kx2 - SQ(omega / hs.cP));
            
            // Intel and GFortran compilers return different branches of the SQRT for negative reals
            // LP: looks like this just means we want the positive branch
            if(kzP.real() == RL(0.0) && kzP.imag() < RL(0.0)) kzP = -kzP;
            f = kzP;
            g = hs.rho;
        }
        
        Refl = -(o.rho * f - J * kz * g) / (o.rho * f + J * kz * g); // complex reflection coef.
        /*
        printf("cS cP rho (%g,%g) (%g,%g) %g\n", hs.cS.real(), hs.cS.imag(),
            hs.cP.real(), hs.cP.imag(), hs.rho);
        printf("kx kz f g Refl (%g,%g) (%g,%g) (%g,%g) (%g,%g) (%g,%g)\n",
            kx.real(), kx.imag(), kz.real(), kz.imag(), f.real(), f.imag(),
            g.real(), g.imag(), Refl.real(), Refl.imag());
        */
        
        if(STD::abs(Refl) < FL(1.0e-5)){ // kill a ray that has lost its energy in reflection
            newPoint.Amp   = FL(0.0);
            newPoint.Phase = oldPoint.Phase;
        }else{
            newPoint.Amp   = STD::abs(Refl) * oldPoint.Amp;
            newPoint.Phase = oldPoint.Phase + STD::atan2(Refl.imag(), Refl.real());
            
            // compute beam-displacement Tindle, Eq. (14)
            // needs a correction to beam-width as well ...
            // LP: most of these variables don't exist, likely very old code
            // if(kz2Sq.real() < FL(0.0)){
            //     rhoW = FL(1.0); // density of water
            //     rhoWSq = rhoW * rhoW;
            //     rhoHSSq = rhoHS * rhoHS;
            //     delta = FL(2.0) * gk * rhoW * rhoS * (kz1Sq - kz2Sq) /
            //         (kz1 * i * kz2 *
            //             (-rhoWSq * kz2Sq + rhoHSSq * kz1Sq));
            //     rv[is+1] = rv[is+1] + delta;
            // }
            
            if(Beam->Type[3] == 'S'){ // beam displacement & width change (Seongil's version)
                ch = oldPoint.c / STD::conj(hs.cP);
                co = oldPoint.t.x * oldPoint.c;
                si = oldPoint.t.y * oldPoint.c;
                ck = omega / oldPoint.c;
                
                a   = FL(2.0) * hs.rho * (FL(1.0) - SQ(ch));
                b   = SQ(co) - SQ(ch);
                d   = SQ(hs.rho) * SQ(si) + b;
                sb  = STD::sqrt(b);
                cco = SQ(co);
                ssi = SQ(si);
                
                if(si != FL(0.0)){
                    delta = a * co / si / (ck * sb * d); // Do we need an abs() on this???
                }else{
                    delta = FL(0.0);
                }
                
                pdelta = delta.real() / (oldPoint.c / co);
                // LP: The spacing in the original version of this formula,
                // the fact that several terms could be factored out to reduce
                // computation, and the repeated divisons, lead me to believe
                // that it may not be correct.
                // Here is the original version with the weird spacing:
                // ddelta = -a / (ck*sb*d) - a*cco / ssi / (ck*sb*d) + a*cco / (ck*b*sb*d)
                //     -a*co / si / (ck*sb*d*d) * (FL(2.0)* SQ(hs.rho) *si*co-FL(2.0)*co*si);
                // Here is a version with things factored better:
                cpx cksbd = ck * sb * d;
                ddelta = a * (cco / (cksbd * b)
                    - (RL(1.0) + (cco / ssi)) / cksbd
                    - FL(2.0) * SQ(co) * (SQ(hs.rho) - RL(1.0)) / (cksbd * d) );
                rddelta = -ddelta.real();
                sddelta = rddelta / STD::abs(rddelta);
                
                // next 3 lines have an update by Diana McCammon to allow a sloping bottom
                // I think the formulas are good, but this won't be reliable because it doesn't have the logic
                // that tracks crossing into new segments after the ray displacement.
                
                theta_bot = STD::atan(tBdry.y / tBdry.x); // bottom angle
                newPoint.x.x = newPoint.x.x + delta.real() * STD::cos(theta_bot); // range displacement
                newPoint.x.y = newPoint.x.y + delta.real() * STD::sin(theta_bot); // depth displacement
                newPoint.tau = newPoint.tau + pdelta; // phase change
                newPoint.q   = newPoint.q + sddelta * rddelta * si * o.ccpx.real() * oldPoint.p; // beam-width change
            }
        }
    }else{
        printf("Reflect2D: Unknown boundary condition type\n");
        bail();
    }
    
    // printf("Reflection amp changed from to %g %g\n", oldPoint.Amp, newPoint.Amp);
}


}
