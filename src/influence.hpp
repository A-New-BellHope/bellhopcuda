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
 * Dalpha, Dbeta: angular spacing between rays
 * freq: source frequency
 * c: nominal sound speed
 * u [LP: 3D: P]: Pressure field (LP: [NRz][Nr])
 * 
 * [LP: 3D only:] mbp: this routine should be eliminated
 * LP: The conversion of intensity to pressure can't be eliminated (moved into
 * the Influence* functions) because it must occur after the summing of
 * contributions from different rays/beams.
 */
template<bool R3D> HOST_DEVICE inline void ScalePressure(
    real Dalpha, real Dbeta, real c, cpx epsilon1, cpx epsilon2, float *r, 
    cpxf *u, int32_t Ntheta, int32_t NRz, int32_t Nr, const char (&RunType)[7], real freq)
{
    real cnst;
    cpx cnst3d;
    
    // Compute scale factor for field
    if constexpr(R3D){
        if(RunType[1] == 'C'){
            // Cerveny Gaussian beams in Cartesian coordinates
            // epsilon is normally imaginary here, so cnst is complex
            real sqrtc = STD::sqrt(c);
            real sqrtc3 = CUBE(sqrtc);
            // put this factor into the beam instead?
            cnst3d = STD::sqrt(epsilon1 * epsilon2) * freq * Dbeta * Dalpha / sqrtc3;
            // LP: This is applied before the intensity to pressure conversion.
            for(int32_t itheta=0; itheta<Ntheta; ++itheta){
                for(int32_t irz=0; irz<NRz; ++irz){
                    for(int32_t ir=0; ir<Nr; ++ir){
                        size_t addr = ((size_t)itheta * NRz + irz) * Nr + ir;
                        u[addr] *= Cpx2Cpxf(cnst3d);
                    }
                }
            }
        }else{
            cnst3d = cpx(FL(1.0), FL(0.0)); // LP: set but not used
        }
    }else{
        IGNORE_UNUSED(Dbeta);
        IGNORE_UNUSED(epsilon1);
        IGNORE_UNUSED(epsilon2);
        if(RunType[1] == 'C' || RunType[1] == 'R'){
            cnst = -Dalpha * STD::sqrt(freq) / c;
        }else{
            cnst = FL(-1.0);
        }
    }
    
    // For incoherent run, convert intensity to pressure
    if(RunType[0] != 'C'){
        for(int32_t itheta=0; itheta<Ntheta; ++itheta){
            for(int32_t irz=0; irz<NRz; ++irz){
                for(int32_t ir=0; ir<Nr; ++ir){
                    size_t addr = ((size_t)itheta * NRz + irz) * Nr + ir;
                    u[addr] = cpxf(STD::sqrt(u[addr].real()), 0.0f);
                }
            }
        }
    }
    
    if constexpr(!R3D){
        // scale and/or incorporate cylindrical spreading
        for(int32_t ir=0; ir<Nr; ++ir){
            real factor;
            if(RunType[3] == 'X'){ // line source
                const float local_pi = 3.14159265f;
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

template<bool R3D> HOST_DEVICE inline void AddToField(
    cpxf *uAllSources, const cpxf &dfield, int32_t itheta, int32_t ir, int32_t iz,
    const InfluenceRayInfo<R3D> &inflray, const Position *Pos)
{
    size_t base = GetFieldAddr(inflray.init.isx, inflray.init.isy, inflray.init.isz,
        itheta, iz, ir, Pos);
    AtomicAddCpx(&uAllSources[base], dfield);
}

template<bool O3D, bool R3D> HOST_DEVICE inline void ApplyContribution(
    cpxf *uAllSources, real cnst, real w, real omega, cpx delay, real phaseInt,
    real RcvrDeclAngle, real RcvrAzimAngle,
    int32_t itheta, int32_t ir, int32_t iz, int32_t is,
    const InfluenceRayInfo<R3D> &inflray, const rayPt<R3D> &point1,
    const Position *Pos, const BeamStructure *Beam,
    EigenInfo *eigen, const ArrInfo *arrinfo)
{
    if constexpr(O3D && !R3D){
        itheta = inflray.init.ibeta;
    }
    if(Beam->RunType[0] == 'E'){
        // eigenrays
        RecordEigenHit(itheta, ir, iz, is, inflray.init, eigen);
    }else if(IsArrivalsRun(Beam)){
        // arrivals
        AddArr<R3D>(itheta, iz, ir, cnst * w, omega, phaseInt, delay, 
            inflray.init, RcvrDeclAngle, RcvrAzimAngle, point1.NumTopBnc, point1.NumBotBnc,
            arrinfo, Pos);
    }else{
        cpxf dfield;
        if(Beam->RunType[0] == 'C'){
            // coherent TL
            dfield = Cpx2Cpxf(cnst * w * STD::exp(-J * (omega * delay - phaseInt)));
            // omega * SQ(n) / (FL(2.0) * SQ(point1.c) * delay)))) // curvature correction [LP: 2D only]
        }else{
            // incoherent/semicoherent TL
            real v = cnst * STD::exp((omega * delay).imag());
            v = SQ(v) * w;
            if(Beam->Type[0] == 'B'){
                // Gaussian beam
                v *= STD::sqrt(FL(2.0) * REAL_PI);
            }
            dfield = cpxf((float)v, 0.0f);
        }
        AddToField<R3D>(uAllSources, dfield, itheta, ir, iz, inflray, Pos);
    }
}

/**
 * Picks the optimum value for epsilon (LP: "beam constant" / "beam initial conditions")
 * 
 * omega: angular frequency
 * c: sound speed
 * gradc: gradient
 * Dangle: angular spacing for ray fan (LP: either Dalpha or Dbeta)
 * rLoop: loop range
 * EpsMultiplier: multiplier to manually adjust result
 * 
 * LP: For 3D, called twice, once for alpha and once for beta
 */
template<bool R3D> HOST_DEVICE inline cpx PickEpsilon(
    char BeamType0, char BeamType1, real omega, real c, vec2 gradc, real angle,
    real Dangle, real rLoop, real EpsMultiplier)
{
    // LP: BUG: Multiple codepaths do not set epsilonOpt, leads to UB
    real halfwidth = R3D ? RL(0.0) : DEBUG_LARGEVAL;
    cpx epsilonOpt = cpx(DEBUG_LARGEVAL, DEBUG_LARGEVAL);
    bool defaultHalfwidth = true, defaultEps = true, zeroEps = false;
    real cz;
    //const char *tag;
    switch(BeamType0){
    case 'C':
    case 'R':
        //tag = R3D ? "Cerveny style beam" : "Paraxial beams";
        switch(BeamType1){
        case 'F':
            //tag = "Space filling beams";
            break;
        case 'M':
            //tag = "Minimum width beams";
            defaultHalfwidth = false;
            halfwidth = STD::sqrt(FL(2.0) * c * FL(1000.0) * rLoop / omega);
            break;
        case 'W':
            //tag = "WKB beams";
            if(R3D){
                printf("Warning, BeamType[1] = W unimplemented in BELLHOP3D\n");
                defaultHalfwidth = defaultEps = false;
                break;
            }
            halfwidth = REAL_MAX;
            cz = gradc.y;
            epsilonOpt = (cz == FL(0.0)) ? RL(1e10) :
                ((-STD::sin(angle) / STD::cos(SQ(angle))) * c * c / cz);
            defaultHalfwidth = defaultEps = false;
            break;
        case 'C':
            if(R3D){
                //tag = "Cerveny style beam";
                printf("Warning, BeamType[1] = C buggy in BELLHOP3D\n");
            }else{
                printf("Warning, BeamType[1] = C not implemented in BELLHOP (2D)\n");
            }
            defaultHalfwidth = defaultEps = false;
            break;
        default:
            printf("Invalid BeamType[1] %c buggily ignored in PickEpsilon\n", BeamType1);
            defaultHalfwidth = defaultEps = false;
        }
        break;
    case 'G':
        if constexpr(R3D){
            //tag = "Geometric beam, hat-shaped, Cart. coord.";
            zeroEps = true;
        }else{
            //tag = "Geometric hat beams";
        }
        break;
    case '^':
        if constexpr(R3D){
            //tag = "Geometric beam, hat-shaped, Cart. coord.";
            zeroEps = true;
        }else{
            printf("Warning, BeamType[1] = ^ not properly handled in BELLHOP (2D)\n");
            defaultHalfwidth = defaultEps = false;
        }
        break;
    case 'g':
        if constexpr(R3D){
            //tag = "Geometric beam, hat-shaped, Ray coord.";
            zeroEps = true;
        }else{
            //tag = "Geometric hat beams";
        }
        break;
    case 'B':
        if constexpr(R3D){
            //tag = "Geometric beam, Gaussian-shaped, Cart. coord.";
            zeroEps = true;
        }else{
            //tag = "Geometric Gaussian beams";
        }
        break;
    case 'b':
        if constexpr(R3D){
            //tag = "Geometric beam, Gaussian-shaped, Ray coord.";
            zeroEps = true;
        }else{
            printf(BHC_PROGRAMNAME ": Geo Gaussian beams in ray-cent. coords. not "
                "implemented in BELLHOP (2D)\n");
            bail();
        }
        break;
    case 'S':
        //tag = "Simple Gaussian beams";
        // LP: Supported here in 3D even though not supported in Influence 3D.
        break;
    default:
        printf("Invalid BeamType[0] %c ignored in PickEpsilon\n", BeamType0);
        defaultHalfwidth = defaultEps = false;
    }
    
    if(zeroEps){
        epsilonOpt = FL(0.0);
    }else if(defaultEps){
        if(defaultHalfwidth){
            halfwidth  = (Dangle == FL(0.0)) ? FL(0.0) : (FL(2.0) / ((omega / c) * Dangle));
        }
        epsilonOpt = J * FL(0.5) * omega * SQ(halfwidth);
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
template<bool R3D> HOST_DEVICE inline bool IsAtCaustic(
    const InfluenceRayInfo<R3D> &inflray, real q, bool qleq0){
    if(qleq0){
        return (q <= RL(0.0) && inflray.qOld >  RL(0.0)) || (q >= RL(0.0) && inflray.qOld <  RL(0.0));
    }else{
        return (q <  RL(0.0) && inflray.qOld >= RL(0.0)) || (q >  RL(0.0) && inflray.qOld <= RL(0.0));
    }
}

/**
 * phase shifts at caustics
 */
template<bool R3D> HOST_DEVICE inline void IncPhaseIfCaustic(
    InfluenceRayInfo<R3D> &inflray, real q, bool qleq0)
{
    if(IsAtCaustic(inflray, q, qleq0)) inflray.phase += REAL_PI / FL(2.0);
}

/**
 * phase shifts at caustics
 * LP: point.phase is discarded if the condition is met, is this correct?
 */
template<bool R3D> HOST_DEVICE inline real FinalPhase(const rayPt<R3D> &point, 
    const InfluenceRayInfo<R3D> &inflray, real q)
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
        printf("Image index %d must be 1, 2, or 3\n", image);
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
HOST_DEVICE inline bool IsDuplicatePoint(const rayPt<false> &point0, const rayPt<false> &point1)
{
    return STD::abs(point1.x.x - point0.x.x) < RL(1.0e3) * spacing(point1.x.x);
}

HOST_DEVICE inline void Compute_eps_pB_qB(cpx &eps, cpx &pB, cpx &qB,
    const rayPt<false> &point, const InfluenceRayInfo<false> &inflray,
    const BeamStructure *Beam)
{
    if(Beam->Type[1] == 'C'){
        eps = J * STD::abs(point.q.x / point.q.y);
    }else{
        eps = inflray.epsilon1;
    }
    pB = point.p.x + eps * point.p.y;
    qB = point.q.x + eps * point.q.y;
}

/**
 * Form gamma
 */
template<bool O3D> HOST_DEVICE inline cpx Compute_gamma(
    const rayPt<false> &point, const cpx &pB, const cpx &qB,
    const Origin<O3D, false> &org, const SSPStructure *ssp, SSPSegState &iSeg)
{
    vec2 rayt = point.c * point.t; // unit tangent
    vec2 rayn = vec2(rayt.y, -rayt.x); // unit normal
    
    SSPOutputs<false> o;
    EvaluateSSP<O3D, false>(point.x, point.t, o, org, ssp, iSeg);
    
    real csq = SQ(o.ccpx.real());
    real cS = glm::dot(o.gradc, rayt);
    real cN = glm::dot(o.gradc, rayn);
    
    real Tr = rayt.x;
    real Tz = rayt.y;
    
    if(qB != RL(0.0)){
        return FL(0.5) * (pB / qB * SQ(Tr) + 
            FL(2.0) * cN / csq * Tz * Tr - cS / csq * SQ(Tz));
    }else{
        return FL(0.0);
    }
}

/**
 * pre-apply some scaling
 */
template<bool R3D> HOST_DEVICE inline rayPt<R3D> ScaleBeamIf3D(
    const rayPt<R3D> &point, const InfluenceRayInfo<R3D> &inflray)
{
    if constexpr(R3D){
        rayPt<true> ret = point;
        // scaling for geometric beams
        ret.q[0][0] *= inflray.rcp_q0;
        ret.q[1][0] *= inflray.rcp_q0;
        ret.q[0][1] *= inflray.rcp_qhat0;
        ret.q[1][0] *= inflray.rcp_qhat0;
    }else{
        return point;
    }
}

template<typename MV> HOST_DEVICE inline real QScalar(const MV &q);
template<> HOST_DEVICE inline real QScalar(const vec2 &q) { return q.x; }
template<> HOST_DEVICE inline real QScalar(const mat2x2 &q) { return glm::determinant(q); }

// 

template<bool O3D, bool R3D> HOST_DEVICE inline void Init_Influence(
    InfluenceRayInfo<R3D> &inflray,
    const rayPt<R3D> &point0, RayInitInfo &rinit, vec2 gradc,
    const Position *Pos, const Origin<O3D, R3D> &org, const SSPStructure *ssp, SSPSegState &iSeg,
    const AnglesStructure *Angles, const FreqInfo *freqinfo, const BeamStructure *Beam)
{
    if constexpr(R3D) IGNORE_UNUSED(ssp);
    
    inflray.init = rinit;
    inflray.freq0 = freqinfo->freq0;
    inflray.omega = FL(2.0) * REAL_PI * inflray.freq0;
    inflray.c0 = point0.c;
    inflray.xs = point0.x;
    // LP: The 5x version is changed to 50x on both codepaths before it is used.
    // inflray.RadMax = FL(5.0) * ccpx.real() / freqinfo->freq0; // 5 wavelength max radius
    inflray.RadMax = FL(50.0) * point0.c / freqinfo->freq0; // 50 wavelength max radius
    inflray.Dalpha = Angles->alpha.d;
    inflray.Dbeta = Angles->beta.d;
    
    const real BeamWindow = RL(4.0); // LP: Integer (!) in 2D
    inflray.BeamWindow = (Beam->Type[0] == 'B' || Beam->Type[0] == 'b') ? BeamWindow : RL(1.0);
    inflray.iBeamWindow2 = SQ(Beam->iBeamWindow);
    
    // LP: Did not have the abs for SGB, but it has been added.
    inflray.Ratio1 = STD::sqrt(STD::abs(STD::cos(rinit.alpha))); // point source
    if constexpr(R3D){
        inflray.Ratio1 *= STD::sqrt(inflray.Dalpha * inflray.Dbeta) / point0.c;
    }else{
        if(Beam->RunType[3] != 'R'){
            inflray.Ratio1 = RL(1.0); // line source
        }
    }
    if(Beam->Type[0] == 'B' || Beam->Type[0] == 'b'){
        // Gaussian beams
        inflray.Ratio1 /= FL(2.0) * REAL_PI; // sqrt(2*pi) represents a sum of Gaussians in free space
    }
    
    if(Beam->Type[0] == 'R' || Beam->Type[0] == 'C'){
        inflray.epsilon1 = PickEpsilon<R3D>(Beam->Type[0], Beam->Type[1], inflray.omega,
            point0.c, gradc, rinit.alpha, Angles->alpha.d, Beam->rLoop, Beam->epsMultiplier);
        if constexpr(R3D){
            inflray.epsilon2 = PickEpsilon<R3D>(Beam->Type[0], Beam->Type[1], inflray.omega,
                point0.c, gradc, rinit.beta, Angles->beta.d, Beam->rLoop, Beam->epsMultiplier);
        }
    }
    
    // LP: For all geometric types
    inflray.rcp_q0 = inflray.Dalpha / point0.c; // Reference for J = q0 / q = q * rcp_q0
    if constexpr(R3D){
        inflray.rcp_qhat0 = STD::abs(STD::cos(inflray.init.alpha)) * inflray.Dbeta / point0.c;
    }
    
    // LP: For all geometric types
    inflray.phase = FL(0.0);
    
    if(Beam->Type[0] == 'S'){
        // LP: For SGB
        inflray.qOld = FL(1.0);
    }else{
        // LP: For Cartesian types
        inflray.qOld = QScalar(point0.q); // used to track KMAH index
    }
    
    // LP: For RayCen types
    if(Beam->Type[0] == 'g'){
        // LP: For hat raycen
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
        inflray.x = VEC23<R3D>(RL(0.0));
        inflray.lastValid = false;
    }else{
        // LP: For Cartesian types
        inflray.x = point0.x;
    }
    
    // LP: For Cerveny
    inflray.kmah = 1;
    
    if(Beam->Type[0] == 'S'){
        // LP: For SGB
        inflray.ir = 0;
    }else if(Beam->Type[0] == 'B' || Beam->Type[0] == 'G'){
        // LP: For Cartesian types
        // what if never satisfied?
        // what if there is a single receiver (ir = -1 possible)
        // LP: ir is always valid, even if it means not meeting the condition.
        real r;
        if constexpr(R3D){
            // LP: Originally something like glm::distance(point0.x[0:1], inflray.xs[0:1])
            // but the ray always starts at the source
            r = RL(0.0);
        }else{
            r = point0.x.x;
        }
        inflray.ir = BinarySearchGT(Pos->Rr, Pos->NRr, 1, 0, r); // find index of first receiver to the right of rA
        if constexpr(!R3D){
            if(point0.t.x < RL(0.0) && inflray.ir > 0) --inflray.ir; // if ray is left-traveling, get the first receiver to the left of rA
        }
    }
    
    if(Beam->Type[0] == 'C'){
        // LP: For Cerveny cart
        if constexpr(R3D){
            printf("Run Type 'C' not supported at this time\n");
            bail();
        }else{
            // LP: Partially supported in Nx2D (O3D but not R3D)
            cpx eps0, pB0, qB0;
            Compute_eps_pB_qB(eps0, pB0, qB0, point0, inflray, Beam);
            inflray.gamma = Compute_gamma<O3D>(point0, pB0, qB0, org, ssp, iSeg);
        }
    }else{
        // LP: not used
        inflray.gamma = RL(0.0);
    }
    
    
}

// LP: Main functions.

/**
 * Paraxial (Cerveny-style) beams in ray-centered coordinates
 */
template<bool O3D> HOST_DEVICE inline bool Step_InfluenceCervenyRayCen(
    const rayPt<false> &point0, const rayPt<false> &point1,
    InfluenceRayInfo<false> &inflray, 
    int32_t is, cpxf *uAllSources,
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
                        printf("Unbounded beam\n");
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
                        contri *= Hermite(n, inflray.RadMax, FL(2.0) * inflray.RadMax);
                        
                        AddToField(uAllSources, Cpx2Cpxf(contri), O3D ? inflray.init.ibeta : 0,
                            ir, iz, inflray, Pos);
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
template<bool O3D> HOST_DEVICE inline bool Step_InfluenceCervenyCart(
    const rayPt<false> &point0, const rayPt<false> &point1,
    InfluenceRayInfo<false> &inflray, 
    int32_t is, cpxf *uAllSources, const BdryType *Bdry, 
    const Origin<O3D, false> &org, const SSPStructure *ssp, SSPSegState &iSeg,
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
    gamma1 = Compute_gamma(point1, pB1, qB1, org, ssp, iSeg);
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
            printf("Unbounded beam\n");
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
                    printf("Invalid Beam->Nimage %d\n", Beam->Nimage);
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
            
            AddToField(uAllSources, Cpx2Cpxf(contri), O3D ? inflray.init.ibeta : 0,
                ir, iz, inflray, Pos);
        }
    }
    
    return true;
}

#if 0
TODO 3D not implemented yet
/**
 * Geometrically-spreading beams with a hat-shaped beam in ray-centered coordinates
 */
template<bool O3D, bool R3D> HOST_DEVICE inline bool Step_InfluenceGeoHatRayCen(
    const rayPt<R3D> &point0, const rayPt<R3D> &point1,
    InfluenceRayInfo<R3D> &inflray,
    int32_t is, cpxf *uAllSources, const Position *Pos, const BeamStructure *Beam,
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
            l = STD::abs(q) * inflray.rcp_q0; // beam radius
            
            if(n < l){ // in beamwindow?
                delay = point0.tau + w * dtau;
                cnst  = scaledAmp / STD::sqrt(STD::abs(q));
                w     = (l - n) / l; // hat function: 1 on center, 0 on edge
                phaseInt = FinalPhase<R3D>(point0, inflray, q);
                
                ApplyContribution<O3D, R3D>(uAllSources,
                    cnst, w, inflray.omega, delay, phaseInt,
                    RcvrDeclAngle, RcvrAzimAngle, itheta, ir, iz, is,
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
#endif

/**
 * Geometric, hat-shaped or Gaussian beams in Cartesian coordintes
 *
 * uAllSources: complex pressure field
 */
template<bool O3D, bool R3D> HOST_DEVICE inline bool Step_InfluenceGeoHatOrGaussianCart(
    bool isGaussian, const rayPt<R3D> &point0, const rayPt<R3D> &point1,
    InfluenceRayInfo<R3D> &inflray, int32_t is, cpxf *uAllSources, const Position *Pos,
    const BeamStructure *Beam, EigenInfo *eigen, const ArrInfo *arrinfo)
{
    // LP: Replaced ScaleBeam in 3D with applying the same scale factors below.
    // This avoids modifying the ray and makes the codepaths more similar.
    
    real rA, rB;
    if constexpr(R3D){
        // LP: 3D updates rA even in the early return conditions below, so
        // it doesn't use inflray.x at all.
        rA = glm::length(vec2(point0.x.x, point0.x.y) - vec2(inflray.xs.x, inflray.xs.y));
        rB = glm::length(vec2(point1.x.x, point1.x.y) - vec2(inflray.xs.x, inflray.xs.y));
        if(STD::abs(rB - rA) <= RL(1e3) * spacing(rA)) return true;
        // LP: This is silly logic, see comments in influence3D.f90
        if(is == 0){
            inflray.ir = (rB > rA) ? 0 : Pos->NRr-1;
        }
    }else{
        // LP: This is different from point0.x.x due to early return for duplicate points.
        rA = inflray.x.x;
        rB = point1.x.x;
    }
    VEC23<R3D> x_ray = point0.x;
    
    // compute normalized tangent (compute it because we need to measure the step length)
    VEC23<R3D> rayt = point1.x - point0.x;
    real rlen = glm::length(rayt);
    // if duplicate point in ray, skip to next step along the ray
    // LP: 2D: and don't update rA (inflray.x) for next time
    real rlen_small = RL(1.0e3) * spacing(point1.x.x);
    if(rlen < rlen_small || (R3D && rlen == rlen_small)) return true;
    rayt /= rlen;
    
    // LP: Ray normals
    VEC23<R3D> rayn1, rayn2; // LP: e1, e2 in 3D; rayn, (none) in 2D
    real RcvrDeclAngle, RcvrAzimAngle, angleXY;
    vec2 n_ray_theta = vec2(-rayt.y, rayt.x); // 3D: normal to the ray in the horizontal receiver plane
    angleXY = RadDeg * STD::atan2(rayt.y, rayt.x);
    if constexpr(R3D){
        RayNormal_unit(rayt, point1.phi, rayn1, rayn2);
        RcvrDeclAngle = RadDeg * STD::atan2(rayt.z, glm::length(n_ray_theta));
        RcvrAzimAngle = angleXY;
    }else{
        rayn1 = n_ray_theta; // unit normal to ray
        IGNORE_UNUSED(rayn2);
        RcvrDeclAngle = angleXY;
        RcvrAzimAngle = DEBUG_LARGEVAL;
    }
    
    // LP: Quantities to be interpolated between steps
    V2M2<R3D> dq   = point1.q   - point0.q;   // LP: dqds in 2D
    cpx       dtau = point1.tau - point0.tau; // LP: dtauds in 2D
    
    // phase shifts at caustics
    real phaseq = QScalar(point0.q);
    IncPhaseIfCaustic(inflray, phaseq, true);
    inflray.qOld = phaseq;
    
    real lambda = isGaussian ? (point0.c / inflray.freq0) : RL(0.0); // local wavelength
    real L_diag; // LP: 3D
    real zmin, zmax; // LP: 2D here, 3D later
    // beam window: kills beams outside exp(RL(-0.5) * SQ(ibwin))
    if constexpr(R3D){
        // beamwidths / LP: Variable values don't carry over to per-receiver beamwidth below
        real l1 = bhc::max(
            glm::length(glm::row(point0.q, 0) * inflray.rcp_q0   ),
            glm::length(glm::row(point1.q, 0) * inflray.rcp_q0   ));
        real l2 = bhc::max(
            glm::length(glm::row(point0.q, 1) * inflray.rcp_qhat0),
            glm::length(glm::row(point1.q, 1) * inflray.rcp_qhat0));
        // worst case is when rectangle is rotated to catch the hypotenuse
        L_diag = STD::sqrt(SQ(l1) + SQ(l2));
    }else{
        real sigma, RadiusMax;
        sigma = bhc::max(STD::abs(point0.q.x), STD::abs(point1.q.x)) * inflray.rcp_q0
            / STD::abs(rayt.x); // beam radius projected onto vertical line
        if(isGaussian){
            // calculate beam width
            sigma     = bhc::max(sigma, bhc::min(FL(0.2) * inflray.freq0 * point1.tau.real(), REAL_PI * lambda));
        }
        RadiusMax = inflray.BeamWindow * sigma; // LP: 1 * sigma for non-Gaussian
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
        IGNORE_UNUSED(L_diag);
    }
    
    // compute beam influence for this segment of the ray
    while(true){
        // is Rr[ir] contained in [rA, rB)? Then compute beam influence
        // LP: Because of the new setup and always incrementing regardless of
        // which direction the ray goes, we only have to check this side.
        if(Pos->Rr[inflray.ir] >= bhc::min(rA, rB) && Pos->Rr[inflray.ir] < bhc::max(rA, rB)){
            
            int32_t itheta = -1;
            do{
                VEC23<R3D> x_rcvr;
                
                ++itheta;
                if constexpr(R3D){
                    // LP: Loop logic
                    if(itheta >= Pos->Ntheta) break;
                    
                    vec2 t_rcvr = Pos->t_rcvr[itheta];
                    SETXY(x_rcvr, XYCOMP(inflray.xs) + (real)Pos->Rr[inflray.ir] * t_rcvr);
                    // normal distance from rcvr to ray segment
                    real m_prime = STD::abs(glm::dot(XYCOMP(x_rcvr) - XYCOMP(x_ray), n_ray_theta));
                    // LP: Commented out in Gaussian
                    if(!isGaussian && m_prime > inflray.BeamWindow * L_diag) continue;
                    
                    // The set of possible receivers is a ring
                    // However, extrapolating the beam backwards produces contributions with s negative and large
                    // We do not want to accept these contributions--- they have the proper range but are 180 degrees
                    // away from this segment of the ray
                    // LP: This s value is not reused below
                    // vector to rcvr dotted into vector to ray point
                    real s = glm::dot(XYCOMP(x_rcvr) - XYCOMP(inflray.xs), XYCOMP(x_ray) - XYCOMP(inflray.xs));
                    if(s < RL(0.0)) continue;
                    
                    // calculate z-limits for the beam (could be pre-cacluated for each itheta)
                    vec2 e_theta = vec2(-t_rcvr.y, t_rcvr.x); // normal to the vertical receiver plane
                    real n_ray_z = rayt.x * e_theta.y - rayt.y * e_theta.x; // normal to the ray in the vertical receiver plane
                    
                    if(STD::abs(n_ray_z) < RL(1e-9)) continue; // avoid divide by zero
                    real L_z = inflray.BeamWindow * L_diag / STD::abs(n_ray_z);
                    
                    zmin = bhc::min(point0.x.z, point1.x.z) - L_z; // min depth of ray segment
                    zmax = bhc::max(point0.x.z, point1.x.z) + L_z; // max depth of ray segment
                }else{
                    IGNORE_UNUSED(itheta);
                    x_rcvr.x = Pos->Rr[inflray.ir];
                }
            
                for(int32_t iz=0; iz<Pos->NRz_per_range; ++iz){
                    int32_t tempiz = iz;
                    if constexpr(!R3D){
                        if(Beam->RunType[4] == 'I') tempiz = inflray.ir; // irregular grid
                    } // else rectilinear grid
                    DEP(x_rcvr) = Pos->Rz[tempiz];
                    if(DEP(x_rcvr) < zmin || DEP(x_rcvr) > zmax) continue;
                    
                    VEC23<R3D> x_rcvr_ray = x_rcvr - x_ray;
                    
                    // linear interpolation of q's
                    real s, n1, n2;
                    s      = glm::dot(x_rcvr_ray, rayt) / rlen; // proportional distance along ray
                    n1     = STD::abs(glm::dot(x_rcvr_ray, rayn1)); // normal distance to ray
                    if constexpr(R3D){
                        n2 = STD::abs(glm::dot(x_rcvr_ray, rayn2)); // normal distance to ray
                    }else{
                        IGNORE_UNUSED(n2);
                    }
                    V2M2<R3D> qInterp = point0.q + s * dq; // interpolated amplitude
                    if constexpr(R3D){
                        qInterp[0][0] *= inflray.rcp_q0;
                        qInterp[1][0] *= inflray.rcp_q0;
                        qInterp[0][1] *= inflray.rcp_qhat0;
                        qInterp[1][1] *= inflray.rcp_qhat0;
                    }
                    real qFinal = QScalar(qInterp); // area of parallelogram formed by ray tube
                    
                    real n1prime, n2prime; // LP: equal to n1 in 2D, rotated & normalized in 3D (were a, b)
                    real sigma, sigma_orig; // beam radius
                    real beamCoordDist;
                    if constexpr(R3D){
                        real l1 = glm::length(glm::row(qInterp, 0));
                        real l2 = glm::length(glm::row(qInterp, 1));
                        if(l1 == FL(0.0) || l2 == FL(0.0)) continue;
                        
                        if(qFinal == FL(0.0)) continue;
                        
                        n1prime = STD::abs((-qInterp[0][1] * n2 + qInterp[1][1] * n1) / qFinal);
                        n2prime = STD::abs(( qInterp[0][0] * n2 - qInterp[1][0] * n1) / qFinal);
                        
                        beamCoordDist = isGaussian ? (n1prime + n2prime) : bhc::max(n1prime, n2prime);
                        sigma = FL(1.0);
                    }else{
                        sigma = sigma_orig = STD::abs(qFinal * inflray.rcp_q0); // LP: called RadiusMax in non-Gaussian
                        if(isGaussian){
                            sigma = bhc::max(sigma, bhc::min(FL(0.2) * inflray.freq0 * point1.tau.real(), REAL_PI * lambda)); // min pi * lambda, unless near
                        }
                        
                        beamCoordDist = n1prime = n1; IGNORE_UNUSED(n2prime);
                    }
                    if(beamCoordDist > inflray.BeamWindow * sigma || 
                        (R3D && beamCoordDist == inflray.BeamWindow * sigma)) continue;
                    
                    cpx delay = point0.tau + s * dtau; // interpolated delay
                    real cnst = inflray.Ratio1 * point1.Amp * STD::sqrt(RL(1.0) / STD::abs(qFinal));
                    real w;
                    if constexpr(R3D){
                        cnst *= point1.c; // LP: From ScaleBeam
                        if(isGaussian){
                            w = STD::exp(FL(-0.5) * (SQ(n1prime) + SQ(n2prime)));
                        }else{
                            w = (FL(1.0) - n1prime) * (FL(1.0) - n2prime);
                        }
                    }else{
                        cnst *= STD::sqrt(point1.c);
                        if(isGaussian){
                            w = STD::exp(FL(-0.5) * SQ(n1prime / sigma)) * (sigma_orig / sigma); // Gaussian decay
                        }else{
                            w = (sigma - n1prime) / sigma; // hat function: 1 on center, 0 on edge
                        }
                    }
                    real phaseInt = FinalPhase<R3D>((!R3D && isGaussian ? point1 : point0),
                        inflray, qFinal);
                    
                    ApplyContribution<O3D, R3D>(uAllSources,
                        cnst, w, inflray.omega, delay, phaseInt,
                        RcvrDeclAngle, RcvrAzimAngle, itheta, inflray.ir, iz, is,
                        inflray, point1, Pos, Beam, eigen, arrinfo);
                }
            } while(R3D);
        }
        
        // receiver not bracketed; bump receiver index, inflray.ir, towards rB
        int32_t irTT;
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
template<bool O3D> HOST_DEVICE inline bool Step_InfluenceSGB(
    const rayPt<false> &point0, const rayPt<false> &point1,
    InfluenceRayInfo<false> &inflray, int32_t is, cpxf *uAllSources,
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
            if(Adeltaz < inflray.RadMax || IsTLRun(Beam)){
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
                ApplyContribution<O3D, false>(uAllSources,
                    cnst, w, inflray.omega, delay, phaseInt,
                    RcvrDeclAngle, RL(0.0), 0, inflray.ir, iz, is,
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
template<bool O3D, bool R3D> HOST_DEVICE inline bool Step_Influence(
    const rayPt<R3D> &point0, const rayPt<R3D> &point1,
    InfluenceRayInfo<R3D> &inflray, int32_t is, cpxf *uAllSources, const BdryType *Bdry,
    const Origin<O3D, R3D> &org, const SSPStructure *ssp, SSPSegState &iSeg,
    const Position *Pos, const BeamStructure *Beam, EigenInfo *eigen, const ArrInfo *arrinfo)
{
    if constexpr(R3D){
        IGNORE_UNUSED(Bdry);
        IGNORE_UNUSED(ssp);
    }
    
    switch(Beam->Type[0]){
    case 'R': 
        if constexpr(R3D){
            printf("Invalid Run Type\n");
            bail();
            return false;
        }else{
            return Step_InfluenceCervenyRayCen<O3D>(
                point0, point1, inflray, is, uAllSources, Bdry, Pos, Beam);
        }
    case 'C':
        if constexpr(O3D){
            printf("Run Type 'C' not supported at this time\n");
            bail();
        }
        if constexpr(R3D){
            // LP: Commented out code to "assemble f, g, h from p-q"
            // using the commented out variables in the 3D ray step
            // return Step_Influence3D(...);
            return false;
        }else{
            // LP: After the error, which does stop the program, Nx2D continues
            // to compute InfluenceCervenyCart! So this has to compile.
            // Nx2D also supports CervenyRayCen above.
            return Step_InfluenceCervenyCart<O3D>(
                point0, point1, inflray, is, uAllSources, Bdry, org, ssp, iSeg, Pos, Beam);
        }
    case 'g':
        /*
        return Step_InfluenceGeoHatRayCen<O3D, R3D>(
            point0, point1, inflray, is, uAllSources, Pos, Beam, eigen, arrinfo);
        */
        printf("3D RayCen not implemented yet\n");
        bail();
        return false;
    case 'S':
        if constexpr(R3D){
            printf("Invalid Run Type\n");
            bail();
            return false;
        }else{
            return Step_InfluenceSGB<O3D>(
                point0, point1, inflray, is, uAllSources, Pos, Beam, eigen, arrinfo);
        }
    case 'B':
    default:
        if constexpr(O3D){
            switch(Beam->Type[0]){
            case 'b':
                printf("GaussianRayCen not implemented yet\n");
                bail();
                return false; /*Step_InfluenceGeoGaussianRayCen(blah);*/
            case 'B':
            case 'G':
            case '^':
                break; //HatOrGaussianCart below
            default:
                printf("Invalid Run Type\n");
                bail();
                return false;
            }
        }
        return Step_InfluenceGeoHatOrGaussianCart<O3D, R3D>(Beam->Type[0] == 'B',
            point0, point1, inflray, is, uAllSources, Pos, Beam, eigen, arrinfo);
    }
}

}
