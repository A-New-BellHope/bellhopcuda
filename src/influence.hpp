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
#include "eigenrays.hpp"
#include "arrivals.hpp"

namespace bhc {

struct InfluenceRayInfo {
    // LP: Variables carried over between iterations.
    real phase, qOld;
    vec2 x;
    bool lastValid;
    real zn, rn;
    int32_t kmah;
    int32_t ir;
    cpx gamma;
    // LP: Constants.
    int32_t isrc, ialpha; // LP: Uniquely identifies the beam
    real q0;
    cpx epsilon; // beam constant
    real freq0, omega;
    real RadMax;
    int32_t iBeamWindow2;
    real SrcDeclAngle; // take-off angle in degrees
    real Dalpha; // angular spacing
    real Ratio1; // scale factor (point source vs. line source)
};

/**
 * Calculates a smoothing function based on the h0 hermite cubic
 * x is the point where the function is to be evaluated
 * returns:
 * [ 0, x1 ] = 1
 * [x1, x2 ] = cubic taper from 1 to 0
 * [x2, inf] = 0
 */
HOST_DEVICE inline real Hermite(real x, real x1, real x2)
{
    real Ax = STD::abs(x);
    if(Ax <= x1){
        return RL(1.0);
    }else if(Ax >= x2){
        return RL(0.0);
    }else{
        real u = (Ax - x1) / (x2 - x1);
        return (RL(1.0) + RL(2.0) * u) * SQ(RL(1.0) - u);
    }
    // ret /= (FL(0.5) * (x1 + x2));
}

/**
 * Scale the pressure field
 * 
 * r: ranges (LP: [Nr])
 * Dalpha: angular spacing between rays
 * freq: source frequency
 * c: nominal sound speed (LP: [NRz][Nr])
 * u: Pressure field
 */
HOST_DEVICE inline void ScalePressure(real Dalpha, real c, float *r, 
    cpxf *u, int32_t NRz, int32_t Nr, const char (&RunType)[7], real freq)
{
    const float local_pi = 3.14159265f;
    
    // Compute scale factor for field
    real cnst;
    if(RunType[1] == 'C' || RunType[1] == 'R'){
        // Cerveny Gaussian beams in Cartesian or Ray-centered coordinates
        cnst = -Dalpha * STD::sqrt(freq) / c;
    }else{
        cnst = FL(-1.0);
    }
    
    // For incoherent run, convert intensity to pressure
    if(RunType[0] != 'C'){
        for(int32_t irz=0; irz<NRz; ++irz){
            for(int32_t ir=0; ir<Nr; ++ir){
                u[irz*Nr+ir] = cpxf((float)STD::sqrt(u[irz*Nr+ir].real()), 0.0f);
            }
        }
    }
    
    // scale and/or incorporate cylindrical spreading
    for(int32_t ir=0; ir<Nr; ++ir){
        real factor;
        if(RunType[3] == 'X'){ // line source
            factor = FL(-4.0) * STD::sqrt(local_pi) * cnst;
        }else{ // point source
            if(r[ir] == 0.0f){
                factor = RL(0.0); // avoid /0 at origin, return pressure = 0
            }else{
                factor = cnst / (real)STD::sqrt(STD::abs(r[ir]));
            }
        }
        for(int32_t irz=0; irz<NRz; ++irz){
            u[irz*Nr+ir] *= (float)factor;
        }
    }
}

/**
 * Checks for a branch cut crossing and updates kmah accordingly
 */
HOST_DEVICE inline void BranchCut(const cpx &q1C, const cpx &q2C,
    const char (&BeamType)[4], int32_t &kmah)
{
    real q1, q2;
    if(BeamType[1] == 'W'){ // WKBeams
        q1 = q1C.real();
        q2 = q2C.real();
    }else{
        if(q2C.real() >= FL(0.0)) return;
        q1 = q1C.imag();
        q2 = q2C.imag();
    }
    if( (q1 < FL(0.0) && q2 >= FL(0.0)) || 
        (q1 > FL(0.0) && q2 <= FL(0.0)) ) kmah = -kmah;
}

HOST_DEVICE inline void ApplyContribution(
    cpxf *u, real cnst, real w, real omega, cpx delay, real phaseInt,
    real SrcDeclAngle, real RcvrDeclAngle, int32_t ir, int32_t iz, int32_t is,
    const InfluenceRayInfo &inflray, const ray2DPt &point1,
    const Position *Pos, const BeamStructure *Beam,
    EigenInfo *eigen, const ArrInfo *arrinfo)
{
    switch(Beam->RunType[0]){
    case 'E':
        // eigenrays
        RecordEigenHit(ir, iz, inflray.isrc, inflray.ialpha, is, eigen);
        break;
    case 'A':
    case 'a':
        // arrivals
        AddArr(omega, inflray.isrc, iz, ir, cnst * w, phaseInt, delay, 
            SrcDeclAngle, RcvrDeclAngle, point1.NumTopBnc, point1.NumBotBnc,
            arrinfo, Pos);
        break;
    case 'C':
        // coherent TL
        AtomicAddCpx(u, Cpx2Cpxf(cnst * w * STD::exp(-J * (omega * delay - phaseInt))));
        // omega * SQ(n) / (FL(2.0) * SQ(point1.c) * delay)))) // curvature correction
        break;
    default:
        // incoherent/semicoherent TL
        real v = cnst * STD::exp((omega * delay).imag());
        v = SQ(v) * w;
        if(Beam->Type[0] == 'B'){
            // Gaussian beam
            v *= STD::sqrt(FL(2.0) * REAL_PI);
        }
        AtomicAddCpx(u, cpxf((float)v, 0.0f));
    }
}

/**
 * Picks the optimum value for epsilon
 * 
 * omega: angular frequency
 * c: sound speed
 * gradc: gradient
 * alpha: angular spacing for ray fan
 * rLoop: loop range
 * EpsMultiplier: multiplier to manually adjust result
 */
HOST_DEVICE inline cpx PickEpsilon(char BeamType0, char BeamType1, real omega, 
    real c, vec2 gradc, real alpha, real Dalpha, real rLoop, real EpsMultiplier)
{
    real halfwidth, cz;
    cpx epsilonOpt;
    //const char *tag;
    switch(BeamType0){
    case 'C':
    case 'R':
        //tag = "Paraxial beams";
        switch(BeamType1){
        case 'F':
            //tag = "Space filling beams";
            halfwidth = FL(2.0) / ((omega / c) * Dalpha);
            epsilonOpt = J * FL(0.5) * omega * SQ(halfwidth);
            break;
        case 'M':
            //tag = "Minimum width beams";
            halfwidth = STD::sqrt(FL(2.0) * c * FL(1000.0) * rLoop / omega);
            epsilonOpt = J * FL(0.5) * omega * SQ(halfwidth);
            break;
        case 'W':
            //tag = "WKB beams";
            halfwidth = REAL_MAX;
            cz = gradc.y;
            if(cz == FL(0.0)){
                epsilonOpt = RL(1e10);
            }else{
                epsilonOpt = (-STD::sin(alpha) / STD::cos(SQ(alpha))) * c * c / cz;
            }
            break;
        default:
            GlobalLog("Invalid BeamType[1]: %c\n", BeamType1);
            bail();
        }
        break;
    case 'G':
    case 'g':
        //tag = "Geometric hat beams";
        halfwidth = FL(2.0) / ((omega / c) * Dalpha);
        epsilonOpt = J * FL(0.5) * omega * SQ(halfwidth);
        break;
    case 'B':
        //tag = "Geometric Gaussian beams";
        halfwidth = FL(2.0) / ((omega / c) * Dalpha);
        epsilonOpt = J * FL(0.5) * omega * SQ(halfwidth);
        break;
    case 'b':
        GlobalLog(BHC_PROGRAMNAME ": Geo Gaussian beams in ray-cent. coords. not "
            "implemented in BELLHOP (and therefore not in " BHC_PROGRAMNAME ")\n");
        bail();
        break;
    case 'S':
        //tag = "Simple Gaussian beams";
        halfwidth = FL(2.0) / ((omega / c) * Dalpha);
        epsilonOpt = J * FL(0.5) * omega * SQ(halfwidth);
        break;
    default:
        GlobalLog("Invalid BeamType[0]: %c\n", BeamType0);
        bail();
    }
    
    /*
    // On first call write info to prt file
    static bool INIFlag = true;
    if(INIFlag){
        PRTFile << "\n" << tag << "\n";
        PRTFile << "halfwidth  = " << halfwidth << "\n";
        PRTFile << "epsilonOpt = " << epsilonOpt << "\n";
        PRTFile << "EpsMult    = " << EpsMultiplier << "\n\n";
        INIFlag = false;
    }
    */
    
    return EpsMultiplier * epsilonOpt;
}

// LP: Helper functions.

/**
 * LP: There are two versions of the phase shift condition used in the BELLHOP
 * code, with the equalities in opposite positions. qleq0 false is only used in
 * SGB.
 */
HOST_DEVICE inline bool IsAtCaustic(const InfluenceRayInfo &inflray, real q, bool qleq0){
    if(qleq0){
        return (q <= RL(0.0) && inflray.qOld >  RL(0.0)) || (q >= RL(0.0) && inflray.qOld <  RL(0.0));
    }else{
        return (q <  RL(0.0) && inflray.qOld >= RL(0.0)) || (q >  RL(0.0) && inflray.qOld <= RL(0.0));
    }
}

/**
 * phase shifts at caustics
 */
HOST_DEVICE inline void IncPhaseIfCaustic(InfluenceRayInfo &inflray, real q, bool qleq0)
{
    if(IsAtCaustic(inflray, q, qleq0)) inflray.phase += REAL_PI / FL(2.0);
}

/**
 * phase shifts at caustics
 * LP: point.phase is discarded if the condition is met, is this correct?
 */
HOST_DEVICE inline real FinalPhase(const ray2DPt &point, 
    const InfluenceRayInfo &inflray, real q)
{
    real phaseInt = point.Phase + inflray.phase;
    if(IsAtCaustic(inflray, q, true)) phaseInt = inflray.phase + REAL_PI / FL(2.0);
    return phaseInt;
}

HOST_DEVICE inline int32_t RToIR(real r, const Position *Pos)
{
    // mbp: should be ", -1);"? [LP: for computation of ir1 / irA]
    return bhc::max(bhc::min((int)((r - Pos->Rr[0]) / Pos->Delta_r), Pos->NRr-1), 0);
}

HOST_DEVICE inline real FlipBeamForImage(real xy, int32_t image, const BdryType *Bdry)
{
    if(image == 1){ // True beam
        return xy;
    }else if(image == 2){ // Surface-reflected beam
        return FL(2.0) * Bdry->Top.hs.Depth - xy;
    }else if(image == 3){ // Bottom-reflected beam
        return FL(2.0) * Bdry->Bot.hs.Depth - xy;
    }else{
        GlobalLog("Image index %d must be 1, 2, or 3\n", image);
        bail();
    }
}

HOST_DEVICE inline void Compute_N_R_IR(real &n, real &r, int32_t &ir,
    real xx, real xy, real zn, real rn, real zR, const Position *Pos)
{
    n = (zR - xy) / zn;
    r = xx + n * rn;
    // following assumes uniform spacing in Pos->r
    ir = RToIR(r, Pos); // index of receiver
}

/**
 * detect and skip duplicate points (happens at boundary reflection)
 */
HOST_DEVICE inline bool IsDuplicatePoint(const ray2DPt &point0, const ray2DPt &point1)
{
    return STD::abs(point1.x.x - point0.x.x) < RL(1.0e3) * spacing(point1.x.x);
}

HOST_DEVICE inline void Compute_eps_pB_qB(cpx &eps, cpx &pB, cpx &qB,
    const ray2DPt &point, const InfluenceRayInfo &inflray, const BeamStructure *Beam)
{
    if(Beam->Type[1] == 'C'){
        eps = J * STD::abs(point.q.x / point.q.y);
    }else{
        eps = inflray.epsilon;
    }
    pB = point.p.x + eps * point.p.y;
    qB = point.q.x + eps * point.q.y;
}

/**
 * Form gamma
 */
HOST_DEVICE inline cpx Compute_gamma(const ray2DPt &point, const cpx &pB, const cpx &qB,
    const InfluenceRayInfo &inflray, const SSPStructure *ssp,
    int32_t &iSegz, int32_t &iSegr)
{
    vec2 rayt = point.c * point.t; // unit tangent
    vec2 rayn = vec2(rayt.y, -rayt.x); // unit normal
    
    cpx ccpx; vec2 gradc; real crr, crz, czz, rho;
    EvaluateSSP(point.x, point.t, ccpx, gradc, crr, crz, czz, rho, inflray.freq0, ssp, iSegz, iSegr);
    
    real csq = SQ(ccpx.real());
    real cS = glm::dot(gradc, rayt);
    real cN = glm::dot(gradc, rayn);
    
    real Tr = rayt.x;
    real Tz = rayt.y;
    
    if(qB != RL(0.0)){
        return FL(0.5) * (pB / qB * SQ(Tr) + 
            FL(2.0) * cN / csq * Tz * Tr - cS / csq * SQ(Tz));
    }else{
        return FL(0.0);
    }
}

// 

HOST_DEVICE inline void Init_Influence(InfluenceRayInfo &inflray,
    const ray2DPt &point0, int32_t isrc, int32_t ialpha, real alpha, vec2 gradc,
    const Position *Pos, const SSPStructure *ssp, int32_t &iSegz, int32_t &iSegr,
    const AnglesStructure *Angles, const FreqInfo *freqinfo, const BeamStructure *Beam)
{
    inflray.isrc = isrc;
    inflray.ialpha = ialpha;
    inflray.freq0 = freqinfo->freq0;
    inflray.omega = FL(2.0) * REAL_PI * inflray.freq0;
    // LP: The 5x version is changed to 50x on both codepaths before it is used.
    // inflray.RadMax = FL(5.0) * ccpx.real() / freqinfo->freq0; // 5 wavelength max radius
    inflray.RadMax = FL(50.0) * point0.c / freqinfo->freq0; // 50 wavelength max radius
    inflray.iBeamWindow2 = SQ(Beam->iBeamWindow);
    inflray.SrcDeclAngle = RadDeg * alpha;
    inflray.Dalpha = Angles->Dalpha;
    if(Beam->RunType[3] == 'R'){
        // LP: Did not have the abs for SGB, but it has been added.
        inflray.Ratio1 = STD::sqrt(STD::abs(STD::cos(alpha))); // point source
    }else{
        inflray.Ratio1 = RL(1.0); // line source
    }
    if(Beam->Type[0] == 'R' || Beam->Type[0] == 'C'){
        inflray.epsilon = PickEpsilon(Beam->Type[0], Beam->Type[1], inflray.omega,
            point0.c, gradc, alpha, Angles->Dalpha, Beam->rLoop, Beam->epsMultiplier);
    }
    
    // LP: For all except Cerveny
    inflray.phase = FL(0.0);
    
    if(Beam->Type[0] == 'S'){
        // LP: For SGB
        inflray.qOld = FL(1.0);
    }else{
        // LP: For GeoHat and Gaussian
        inflray.qOld = point0.q.x; // used to track KMAH index
    }
    
    // LP: For GeoHat and Gaussian
    inflray.q0 = point0.c / inflray.Dalpha; // Reference for J = q0 / q
    
    // LP: For RayCen types
    if(Beam->Type[0] == 'g'){
        // LP: For GeoHat
        inflray.zn = -point0.t.x * point0.c;
        inflray.rn =  point0.t.y * point0.c;
        inflray.x = point0.x;
        inflray.lastValid = STD::abs(inflray.zn) >= RL(1e-6);
    }else if(Beam->Type[0] == 'R'){
        // LP: For Cerveny
        // mbp: This logic means that the first step along the ray is skipped
        // which is a problem if deltas is very large, e.g. isospeed problems
        // I [mbp] fixed this in InfluenceGeoHatRayCen
        inflray.zn = RL(1.0);
        inflray.rn = RL(0.0);
        inflray.x = vec2(RL(0.0), RL(0.0));
        inflray.lastValid = false;
    }else{
        // LP: For geo Cartesian types
        inflray.x = point0.x;
    }
    
    // LP: For Cerveny
    inflray.kmah = 1;
    
    if(Beam->Type[0] == 'S'){
        // LP: For SGB
        inflray.ir = 0;
    }else if(Beam->Type[0] == 'B' || Beam->Type[0] == 'G'){
        // LP: For GeoHat and Gaussian
        // what if never satisfied?
        // what if there is a single receiver (ir = -1 possible)
        // LP: ir is always valid, even if it means not meeting the condition.
        inflray.ir = BinarySearchGT(Pos->Rr, Pos->NRr, 1, 0, point0.x.x); // find index of first receiver to the right of rA
        if(point0.t.x < RL(0.0) && inflray.ir > 0) --inflray.ir; // if ray is left-traveling, get the first receiver to the left of rA
    }
    
    if(Beam->Type[0] == 'C'){
        // LP: For Cerveny cart
        cpx eps0, pB0, qB0;
        Compute_eps_pB_qB(eps0, pB0, qB0, point0, inflray, Beam);
        inflray.gamma = Compute_gamma(point0, pB0, qB0, inflray, ssp, iSegz, iSegr);
    }else{
        // LP: not used
        inflray.gamma = RL(0.0);
    }
}

// LP: Main functions.

/**
 * Paraxial (Cerveny-style) beams in ray-centered coordinates
 */
HOST_DEVICE inline bool Step_InfluenceCervenyRayCen(
    const ray2DPt &point0, const ray2DPt &point1, InfluenceRayInfo &inflray, 
    int32_t is, cpxf *u,
    const BdryType *Bdry, const Position *Pos, const BeamStructure *Beam)
{
    IGNORE_UNUSED(is);
    
    cpx eps0, eps1, pB0, pB1, qB0, qB1, gamma0, gamma1;
    real zn, rn, zR;
    // need to add logic related to NRz_per_range
    
    // During reflection imag(q) is constant and adjacent normals cannot bracket
    // a segment of the TL line, so no special treatment is necessary
    
    Compute_eps_pB_qB(eps0, pB0, qB0, point0, inflray, Beam);
    Compute_eps_pB_qB(eps1, pB1, qB1, point1, inflray, Beam);
    gamma0 = pB0 / qB0;
    gamma1 = pB1 / qB1;
    
    // ray normal based on tangent with c(s) scaling
    zn = -point1.t.x * point1.c;
    rn =  point1.t.y * point1.c;
    // If normal parallel to TL-line, skip to next step on ray
    // LP: Changed from (the FORTRAN equivalent of) REAL_MINPOS as this is the
    // smallest positive floating point number, which would be equivalent to
    // just zn == 0.0.
    if(STD::abs(zn) < REAL_EPSILON) return true;
    
    // detect and skip duplicate points (happens at boundary reflection)
    if(IsDuplicatePoint(point0, point1)){
        inflray.lastValid = true;
        inflray.x = point1.x;
        inflray.zn = zn;
        inflray.rn = rn;
        return true;
    }
    
    // compute KMAH index
    // Following is incorrect for 'Cerveny'-style beamwidth (narrow as possible)
    int32_t old_kmah = inflray.kmah;
    BranchCut(qB0, qB1, Beam->Type, inflray.kmah);
    
    for(int32_t iz=0; iz<Pos->NRz_per_range; ++iz){
        zR = Pos->Rz[iz];
        
        for(int32_t image=1; image <= Beam->Nimage; ++image){
            real Polarity = (image == 2) ? RL(-1.0) : RL(1.0);
            real nA, rA, nB, rB;
            int32_t ir1, ir2; // LP: mbp switches from A/B naming to 1/2 here.
            Compute_N_R_IR(nB, rB, ir2, point1.x.x, FlipBeamForImage(point1.x.y, image, Bdry),
                zn, rn * Polarity, zR, Pos);
            Compute_N_R_IR(nA, rA, ir1, inflray.x.x, FlipBeamForImage(inflray.x.y, image, Bdry),
                inflray.zn, inflray.rn * Polarity, zR, Pos);
            
            if(inflray.lastValid && ir1 < ir2){
                for(int32_t ir=ir1+1; ir<=ir2; ++ir){
                    real w, n, nSq, c;
                    cpx q, gamma, tau, contri;
                    w     = (Pos->Rr[ir] - rA) / (rB - rA);
                    q     = qB0    + w * (qB1 - qB0);
                    gamma = gamma0 + w * (gamma1 - gamma0);
                    n     = nA     + w * (nB - nA);
                    nSq   = SQ(n);
                    if(gamma.imag() > 0){
                        #ifndef BHC_USE_FLOATS
                        GlobalLog("Unbounded beam\n");
                        #endif
                        continue;
                    }
                    
                    if(FL(-0.5) * inflray.omega * gamma.imag() * nSq < inflray.iBeamWindow2){ // Within beam window?
                        c   = point0.c;
                        tau = point0.tau + w * (point1.tau - point0.tau);
                        contri = inflray.Ratio1 * point1.Amp * STD::sqrt(c * STD::abs(eps1) / q) *
                            STD::exp(-J * (inflray.omega * (tau + FL(0.5) * gamma * nSq) - point1.Phase));
                        
                        cpx P_n = -J * inflray.omega * gamma * n * contri;
                        cpx P_s = -J * inflray.omega / c         * contri;
                        switch(Beam->Component){
                        case 'P': // pressure
                            break;
                        case 'V': // vertical component
                            contri = c * (P_n * point1.t.x + P_s * point1.t.y);
                            break;
                        case 'H': // horizontal component
                            contri = c * (-P_n * point1.t.y + P_s * point1.t.x);
                            break;
                        }
                        
                        int32_t kmah = old_kmah;
                        BranchCut(qB0, q, Beam->Type, kmah); // Get correct branch of STD::sqrt
                        
                        if(kmah < 0)   contri = -contri;
                        if(image == 2) contri = -contri;
                        
                        if(Beam->RunType[0] == 'I' || Beam->RunType[0] == 'S'){ // Incoherent or Semi-coherent TL
                            contri = contri * STD::conj(contri);
                        }
                        
                        AtomicAddCpx(&u[iz*Pos->NRr+ir], Cpx2Cpxf(
                            Hermite(n, inflray.RadMax, FL(2.0) * inflray.RadMax) * contri));
                    }
                }
            }
        }
    }
    
    inflray.lastValid = true;
    inflray.x = point1.x;
    inflray.zn = zn;
    inflray.rn = rn;
    return true;
}

/**
 * Paraxial (Cerveny-style) beams in Cartesian coordinates
 */
HOST_DEVICE inline bool Step_InfluenceCervenyCart(
    const ray2DPt &point0, const ray2DPt &point1, InfluenceRayInfo &inflray, 
    int32_t is, cpxf *u, 
    const BdryType *Bdry, const SSPStructure *ssp, int32_t &iSegz, int32_t &iSegr, 
    const Position *Pos, const BeamStructure *Beam)
{
    cpx eps0, eps1, pB0, pB1, qB0, qB1, gamma0, gamma1;
    real zR;
    // need to add logic related to NRz_per_range
    
    // During reflection imag(q) is constant and adjacent normals cannot bracket
    // a segment of the TL line, so no special treatment is necessary
    
    Compute_eps_pB_qB(eps0, pB0, qB0, point0, inflray, Beam);
    Compute_eps_pB_qB(eps1, pB1, qB1, point1, inflray, Beam);
    
    // Form gamma and KMAH index
    // Treatment of KMAH index is incorrect for 'Cerveny' style beam width BeamType
    
    gamma0 = inflray.gamma;
    gamma1 = Compute_gamma(point1, pB1, qB1, inflray, ssp, iSegz, iSegr);
    inflray.gamma = gamma1;
        
    int32_t old_kmah = inflray.kmah;
    BranchCut(qB0, qB1, Beam->Type, inflray.kmah);
    
    if(is == 0) return true; // LP: Skips the first valid pair.
    // LP: Assumes rays may never travel left.
    if(point1.x.x > Pos->Rr[Pos->NRr-1]){
        return false; // LP: Terminates ray.
    }
    real rA = point0.x.x;
    real rB = point1.x.x;
    if(IsDuplicatePoint(point0, point1)) return true; // don't process duplicate points
    
    // Compute upper index on rcvr line
    // Assumes r is a vector of equally spaced points
    int32_t irA = RToIR(rA, Pos);
    int32_t irB = RToIR(rB, Pos);
    
    if(irA >= irB) return true;
    
    for(int32_t ir=irA+1; ir<=irB; ++ir){
        real w, c;
        vec2 x, rayt;
        cpx q, tau, gamma, cnst;                
        w     = (Pos->Rr[ir] - rA) / (rB - rA);
        x     = point0.x   + w * (point1.x   - point0.x);
        rayt  = point0.t   + w * (point1.t   - point0.t);
        c     = point0.c   + w * (point1.c   - point0.c);
        q     = qB0        + w * (qB1        - qB0);
        tau   = point0.tau + w * (point1.tau - point0.tau);
        gamma = gamma0     + w * (gamma1     - gamma0);
        
        if(gamma.imag() > FL(0.0)){
            #ifndef BHC_USE_FLOATS
            GlobalLog("Unbounded beam\n");
            #endif
            continue;
        }
        
        cnst = inflray.Ratio1 * STD::sqrt(c * STD::abs(eps0) / q);
        
        // Get correct branch of STD::sqrt
        int32_t kmah = old_kmah;
        BranchCut(qB0, q, Beam->Type, kmah);
        if(kmah < 0) cnst = -cnst;
        
        for(int32_t iz=0; iz<Pos->NRz_per_range; ++iz){
            zR = Pos->Rz[iz];
            
            cpx contri = FL(0.0);
            for(int32_t image=1; image <= Beam->Nimage; ++image){
                real deltaz, Polarity;
                if(image == 1){ // True beam
                    deltaz = zR - x.y;
                    Polarity = RL(1.0);
                }else if(image == 2){ // Surface reflected beam
                    deltaz = -zR + FL(2.0) * Bdry->Top.hs.Depth - x.y;
                    Polarity = RL(-1.0);
                }else if(image == 3){ // Bottom  reflected beam
                    deltaz = -zR + FL(2.0) * Bdry->Bot.hs.Depth - x.y;
                    Polarity = RL(1.0); // assumes rigid bottom
                }else{
                    GlobalLog("Invalid Beam->Nimage %d\n", Beam->Nimage);
                    bail();
                }
                if(inflray.omega * gamma.imag() * SQ(deltaz) < inflray.iBeamWindow2){
                    contri += Polarity * point1.Amp * 
                        Hermite(deltaz, inflray.RadMax, FL(2.0) * inflray.RadMax) *
                        STD::exp(-J * (inflray.omega * (
                            tau + rayt.y * deltaz + gamma * SQ(deltaz)) - point1.Phase));
                }
            }
            
            if(Beam->RunType[0] == 'C'){ // coherent
                contri = cnst * contri;
            }else if(Beam->RunType[0] == 'I' || Beam->RunType[0] == 'S'){
                contri = cnst * contri;
                contri = contri * STD::conj(contri);
            }
            AtomicAddCpx(&u[iz*Pos->NRr+ir], Cpx2Cpxf(contri));
        }
    }
    
    return true;
}

/**
 * Geometrically-spreading beams with a hat-shaped beam in ray-centered coordinates
 */
HOST_DEVICE inline bool Step_InfluenceGeoHatRayCen(
    const ray2DPt &point0, const ray2DPt &point1, InfluenceRayInfo &inflray,
    int32_t is, cpxf *u, const Position *Pos, const BeamStructure *Beam,
    EigenInfo *eigen, const ArrInfo *arrinfo)
{
    real RcvrDeclAngle;
    real dq, q, w, n, l, cnst, phaseInt;
    real zn, rn, nA, nB, rA, rB;
    int32_t irA, irB;
    cpx dtau, delay;
    
    dq = point1.q.x - point0.q.x;
    dtau = point1.tau - point0.tau;
    
    // ray normal based on tangent with c(s) scaling
    zn = -point1.t.x * point1.c;
    rn =  point1.t.y * point1.c;
    if(STD::abs(zn) < RL(1e-10)) return true;
    
    if(IsDuplicatePoint(point0, point1)){
        inflray.lastValid = true;
        inflray.x = point1.x;
        inflray.zn = zn;
        inflray.rn = rn;
        return true;
    }
    
    RcvrDeclAngle = RadDeg * STD::atan2(point1.t.y, point1.t.x);
    
    // During reflection imag(q) is constant and adjacent normals cannot bracket
    // a segment of the TL line, so no special treatment is necessary
    
    real scaledAmp = inflray.Ratio1 * STD::sqrt(point1.c) * point1.Amp;
    
    for(int32_t iz=0; iz<Pos->NRz_per_range; ++iz){
        real zR = Pos->Rz[iz];
        
        if(!inflray.lastValid){
            nA = RL(1e10);
            rA = RL(1e10);
            irA = 0;
        }else{
            Compute_N_R_IR(nA, rA, irA, inflray.x.x, inflray.x.y, 
                inflray.zn, inflray.rn, zR, Pos);
        }
        Compute_N_R_IR(nB, rB, irB, point1.x.x, point1.x.y, zn, rn, zR, Pos);
        
        if(irA == irB) continue;
        
        q = point0.q.x;
        IncPhaseIfCaustic(inflray, q, true);
        inflray.qOld = q;
        
        // *** Compute contributions to bracketed receivers ***
        
        for(int32_t ir = bhc::min(irA, irB) + 1; ir <= bhc::max(irA, irB); ++ir){
            w = (Pos->Rr[ir] - rA) / (rB - rA);
            n = STD::abs(nA + w * (nB - nA));
            q = point0.q.x + w * dq; // interpolated amplitude
            l = STD::abs(q) / inflray.q0; // beam radius
            
            if(n < l){ // in beamwindow?
                delay = point0.tau + w * dtau;
                cnst  = scaledAmp / STD::sqrt(STD::abs(q));
                w     = (l - n) / l; // hat function: 1 on center, 0 on edge
                phaseInt = FinalPhase(point0, inflray, q);
                
                ApplyContribution(&u[iz*Pos->NRr+ir], 
                    cnst, w, inflray.omega, delay, phaseInt,
                    inflray.SrcDeclAngle, RcvrDeclAngle, ir, iz, is,
                    inflray, point1, Pos, Beam, eigen, arrinfo);
            }
        }
    }
    
    inflray.lastValid = true;
    inflray.x = point1.x;
    inflray.zn = zn;
    inflray.rn = rn;
    return true;
}

/**
 * Geometric, hat-shaped or Gaussian beams in Cartesian coordintes
 *
 * u: complex pressure field
 */
HOST_DEVICE inline bool Step_InfluenceGeoHatOrGaussianCart(
    bool isGaussian, 
    const ray2DPt &point0, const ray2DPt &point1, InfluenceRayInfo &inflray, 
    int32_t is, cpxf *u, const Position *Pos, const BeamStructure *Beam,
    EigenInfo *eigen, const ArrInfo *arrinfo)
{
    const int32_t BeamWindow = 4; // beam window: kills beams outside e**(-0.5 * ibwin**2 )
    
    vec2 x_ray, rayt, rayn, x_rcvr;
    real rA, rB, RcvrDeclAngle;
    real rlen, q, s, n, cnst, w, phaseInt, RadiusMax, zmin, zmax, sigma, lambda, a, dqds;
    cpx dtauds, delay;
    int32_t irTT;
    
    rA = inflray.x.x; // LP: This is different from point0.x.x due to early return for duplicate points.
    rB = point1.x.x;
    x_ray = point0.x;
    
    // compute normalized tangent (compute it because we need to measure the step length)
    rayt = point1.x - point0.x;
    rlen = glm::length(rayt);
    // if duplicate point in ray, skip to next step along the ray
    // LP: and don't update rA (inflray.x) for next time
    if(rlen < RL(1.0e3) * spacing(point1.x.x)) return true;
    rayt /= rlen;
    rayn = vec2(-rayt.y, rayt.x); // unit normal to ray
    RcvrDeclAngle = RadDeg * STD::atan2(rayt.y, rayt.x);
    
    dqds   = point1.q.x - point0.q.x;
    dtauds = point1.tau - point0.tau;
    
    q = point0.q.x;
    IncPhaseIfCaustic(inflray, q, true);
    inflray.qOld = q;
    
    sigma = bhc::max(STD::abs(point0.q.x), STD::abs(point1.q.x)) / 
        (inflray.q0 * STD::abs(rayt.x)); // beam radius projected onto vertical line
    if(isGaussian){
        // calculate beam width
        lambda    = point0.c / inflray.freq0;
        sigma     = bhc::max(sigma, bhc::min(FL(0.2) * inflray.freq0 * point1.tau.real(), REAL_PI * lambda));
        RadiusMax = BeamWindow * sigma;
    }else{
        lambda    = RL(0.0); // LP: compiler incorrectly complains maybe uninitialized
        RadiusMax = sigma;
    }
    
    // depth limits of beam
    // LP: For rays shot at exactly 60 degrees, they will hit this edge case.
    // This is a sharp edge--the handling on each side of this edge may be
    // significantly different. So, moved the edge away from the round number.
    if(STD::abs(rayt.x) > FL(0.50001)){ // shallow angle ray
        zmin = bhc::min(point0.x.y, point1.x.y) - RadiusMax;
        zmax = bhc::max(point0.x.y, point1.x.y) + RadiusMax;
    }else{ // steep angle ray
        zmin = -REAL_MAX;
        zmax =  REAL_MAX;
    }
    
    // compute beam influence for this segment of the ray
    while(true){
        // is Rr[ir] contained in [rA, rB)? Then compute beam influence
        // LP: Because of the new setup and always incrementing regardless of
        // which direction the ray goes, we only have to check this side.
        if(Pos->Rr[inflray.ir] >= bhc::min(rA, rB) && Pos->Rr[inflray.ir] < bhc::max(rA, rB)){
            
            for(int32_t iz=0; iz<Pos->NRz_per_range; ++iz){
                if(Beam->RunType[4] == 'I'){
                    x_rcvr = vec2(Pos->Rr[inflray.ir], Pos->Rz[inflray.ir]); // irregular grid
                }else{
                    x_rcvr = vec2(Pos->Rr[inflray.ir], Pos->Rz[iz]); // rectilinear grid
                }
                if(x_rcvr.y < zmin || x_rcvr.y > zmax) continue;
                
                s     = glm::dot(x_rcvr - x_ray, rayt) / rlen; // proportional distance along ray
                n     = STD::abs(glm::dot(x_rcvr - x_ray, rayn)); // normal distance to ray
                q     = point0.q.x + s * dqds; // interpolated amplitude
                sigma = STD::abs(q / inflray.q0); // beam radius
                real beamWCompare;
                if(isGaussian){
                    sigma = bhc::max(sigma, bhc::min(FL(0.2) * inflray.freq0 * point1.tau.real(), REAL_PI * lambda)); // min pi * lambda, unless near
                    beamWCompare = BeamWindow * sigma;
                }else{
                    RadiusMax = sigma;
                    beamWCompare = RadiusMax;
                }
                
                if(n < beamWCompare){ // Within beam window?
                    delay = point0.tau + s * dtauds; // interpolated delay
                    if(STD::abs(q) == RL(0.0)){
                        vec2 dx = x_rcvr - x_ray;
                        GlobalLog("Divide by q=zero, point0.q.x %f, s %f, dqds %f, iz %d, ir %d, is %d, rayt (%f,%f), rlen %f, x_rcvr - x_ray (%f,%f), x_ray (%f,%f)\n",
                            point0.q.x, s, dqds, iz, inflray.ir, is, rayt.x, rayt.y, rlen, dx.x, dx.y, x_ray.x, x_ray.y);
                        bail();
                    }
                    cnst  = inflray.Ratio1 * STD::sqrt(point1.c / STD::abs(q)) * point1.Amp;
                    if(isGaussian){
                        // sqrt( 2 * pi ) represents a sum of Gaussians in free space
                        cnst /= STD::sqrt(FL(2.0) * REAL_PI);
                        a = STD::abs(inflray.q0 / q);
                        w = STD::exp(FL(-0.5) * SQ(n / sigma)) / (sigma * a); // Gaussian decay
                    }else{
                        w = (RadiusMax - n) / RadiusMax; // hat function: 1 on center, 0 on edge
                    }
                    phaseInt = FinalPhase((isGaussian ? point1 : point0), inflray, q);
                    
                    ApplyContribution(&u[iz*Pos->NRr+inflray.ir],
                        cnst, w, inflray.omega, delay, phaseInt,
                        inflray.SrcDeclAngle, RcvrDeclAngle, inflray.ir, iz, is,
                        inflray, point1, Pos, Beam, eigen, arrinfo);
                }
            }
        }
        
        // receiver not bracketed; bump receiver index, inflray.ir, towards rB
        if(rB > Pos->Rr[inflray.ir]){
            if(inflray.ir >= Pos->NRr - 1) break; // go to next step on ray
            irTT = inflray.ir + 1; // bump right
            if(Pos->Rr[irTT] >= rB) break; // go to next step on ray
        }else{
            if(inflray.ir <= 0) break; // go to next step on ray
            irTT = inflray.ir - 1; // bump left
            if(Pos->Rr[irTT] <= rB) break; // go to next step on ray
        }
        inflray.ir = irTT;
    }
    
    inflray.x = point1.x;
    return true;
}

/**
 * Bucker's Simple Gaussian Beams in Cartesian coordinates
 */
HOST_DEVICE inline bool Step_InfluenceSGB(
    const ray2DPt &point0, const ray2DPt &point1, InfluenceRayInfo &inflray,
    int32_t is, cpxf *u,
    const Position *Pos, const BeamStructure *Beam,
    EigenInfo *eigen, const ArrInfo *arrinfo)
{
    real w;
    vec2 x, rayt;
    cpx tau;
    real RcvrDeclAngle = RadDeg * STD::atan2(point1.t.y, point1.t.x);
    
    const real beta = FL(0.98); // Beam Factor
    real a = FL(-4.0) * STD::log(beta) / SQ(inflray.Dalpha);
    real cn = inflray.Dalpha * STD::sqrt(a / REAL_PI);
    real rA = point0.x.x;
    
    real rB = point1.x.x;
    
    real q = point0.q.x;
    IncPhaseIfCaustic(inflray, q, false);
    inflray.qOld = q;
    
    // Loop over bracketed receiver ranges
    // LP: BUG: This way of setting up the loop (which matches the FORTRAN
    // implementation) assumes the ray always travels towards positive R, which
    // is not true for certain bathymetries (or for rays simply shot backwards,
    // which would also crash during the setup, see Init_Influence).
    while(STD::abs(rB - rA) > RL(1.0e3) * spacing(rA) && rB > Pos->Rr[inflray.ir]){
        w    = (Pos->Rr[inflray.ir] - rA) / (rB - rA);
        x    = point0.x   + w * (point1.x   - point0.x);
        rayt = point0.t   + w * (point1.t   - point0.t);
        q    = point0.q.x + w * (point1.q.x - point0.q.x);
        tau  = point0.tau + w * (point1.tau - point0.tau);
        
        // following is incorrect because ray doesn't always use a step of deltas
        // LP: The while ignores extremely small steps, but those small steps
        // still increment is, so the later ray segments still treat it as if
        // all steps leading up to them were of size deltas.
        real sint = ((real)(is + 1) + w) * Beam->deltas;
        
        IncPhaseIfCaustic(inflray, q, false);
        
        for(int32_t iz=0; iz<Pos->NRz_per_range; ++iz){
            real deltaz = Pos->Rz[iz] - x.y; // ray to rcvr distance
            // LP: Reinstated this condition for eigenrays and arrivals, as
            // without it every ray would be an eigenray / arrival.
            real Adeltaz = STD::abs(deltaz);
            if(Adeltaz < inflray.RadMax || Beam->RunType[0] == 'C'
                    || Beam->RunType[0] == 'I' || Beam->RunType[0] == 'S'){
                // LP: Changed to use ApplyContribution in order to support 
                // incoherent, semi-coherent, and arrivals.
                real cpa = STD::abs(deltaz * (rB - rA)) / STD::sqrt(SQ(rB - rA) 
                    + SQ(point1.x.y - point0.x.y));
                real ds = STD::sqrt(SQ(deltaz) - SQ(cpa));
                real sx1 = sint + ds;
                real thet = STD::atan(cpa / sx1);
                cpx delay = tau + rayt.y * deltaz;
                real cnst = inflray.Ratio1 * cn * point1.Amp / STD::sqrt(sx1);
                w = STD::exp(-a * SQ(thet));
                real phaseInt = point1.Phase + inflray.phase;
                ApplyContribution(&u[iz*Pos->NRr + inflray.ir],
                    cnst, w, inflray.omega, delay, phaseInt,
                    inflray.SrcDeclAngle, RcvrDeclAngle, inflray.ir, iz, is,
                    inflray, point1, Pos, Beam, eigen, arrinfo);
            }
        }
        
        inflray.qOld = q;
        ++inflray.ir;
        if(inflray.ir >= Pos->NRr) return false;
    }
    
    return true;
}

/**
 * LP: Returns whether to continue the ray.
 */
HOST_DEVICE inline bool Step_Influence(
    const ray2DPt &point0, const ray2DPt &point1, InfluenceRayInfo &inflray, 
    int32_t is, cpxf *u,
    const BdryType *Bdry, const SSPStructure *ssp, int32_t &iSegz, int32_t &iSegr, 
    const Position *Pos, const BeamStructure *Beam, EigenInfo *eigen, const ArrInfo *arrinfo)
{
    switch(Beam->Type[0]){
    case 'R': return Step_InfluenceCervenyRayCen(
            point0, point1, inflray, is, u, Bdry, Pos, Beam);
    case 'C': return Step_InfluenceCervenyCart(
            point0, point1, inflray, is, u, Bdry, ssp, iSegz, iSegr, Pos, Beam);
    case 'g': return Step_InfluenceGeoHatRayCen(
            point0, point1, inflray, is, u, Pos, Beam, eigen, arrinfo);
    case 'S': return Step_InfluenceSGB(
            point0, point1, inflray, is, u, Pos, Beam, eigen, arrinfo);
    case 'B':
    default:  return Step_InfluenceGeoHatOrGaussianCart(Beam->Type[0] == 'B',
            point0, point1, inflray, is, u, Pos, Beam, eigen, arrinfo);
    }
}

}
