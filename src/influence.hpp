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
#include "runtype.hpp"
#include "eigenrays.hpp"
#include "arrivals.hpp"

// #define INFL_DEBUGGING_ITHETA 0
// #define INFL_DEBUGGING_IZ 52
// #define INFL_DEBUGGING_IR 230

namespace bhc {

////////////////////////////////////////////////////////////////////////////////
// General helper functions
////////////////////////////////////////////////////////////////////////////////

/**
 * detect and skip duplicate points (happens at boundary reflection)
 */
template<bool R3D> HOST_DEVICE inline bool IsSmallValue(
    real delta, real ref, bool const_thresh)
{
    const real threshold = const_thresh ? FL(1.0e-4) : RL(1.0e3) * spacing(ref);
    if constexpr(R3D)
        return delta <= threshold;
    else
        return delta < threshold;
}

template<bool R3D> HOST_DEVICE inline bool IsDuplicatePoint(
    real a, real b, bool const_thresh)
{
    return IsSmallValue<R3D>(STD::abs(b - a), b, const_thresh);
}

HOST_DEVICE inline bool IsDuplicatePoint(
    const rayPt<false> &point0, const rayPt<false> &point1, bool const_thresh)
{
    return IsDuplicatePoint<false>(point0.x.x, point1.x.x, const_thresh);
}

HOST_DEVICE inline bool IsDuplicatePoint(
    const rayPt<true> &point0, const rayPt<true> &point1, bool const_thresh)
{
    return IsSmallValue<true>(glm::length(point1.x - point0.x), point1.x.x, const_thresh);
}

////////////////////////////////////////////////////////////////////////////////
// Geometrical-only helper functions
////////////////////////////////////////////////////////////////////////////////

/**
 * index of nearest rcvr before normal
 * Compute upper index on rcvr line
 * Assumes Pos->Rr is a vector of equally spaced points
 */
HOST_DEVICE inline int32_t RToIR(real r, const Position *Pos)
{
    real temp = (r - Pos->Rr[0]) / Pos->Delta_r;
    // LP: Added snapping to deal with floating-point error at the int boundaries.
    if(STD::abs(temp - STD::round(temp)) < RL(1e-6)) { temp = STD::round(temp); }
    // mbp: should be ", -1);"? [LP: for 2D computation of ir1 / irA]
    return bhc::max(bhc::min((int)temp, Pos->NRr - 1), 0);
}

template<typename CFG, bool O3D, bool R3D> HOST_DEVICE inline void AdjustSigma(
    real &sigma, const rayPt<R3D> &point0, const rayPt<R3D> &point1,
    const InfluenceRayInfo<R3D> &inflray)
{
    if constexpr(!CFG::infl::IsGaussianGeom()) return;
    real lambda = point0.c / inflray.freq0; // local wavelength
    // min pi * lambda, unless near
    sigma = bhc::max(
        sigma, bhc::min(FL(0.2) * inflray.freq0 * point1.tau.real(), REAL_PI * lambda));
}

/**
 * mbp: sqrt(2*pi) represents a sum of Gaussians in free space
 * LP: Can't be constexpr as sqrt is not without GCC extensions
 */
template<bool R3D> HOST_DEVICE inline real GaussScaleFactor()
{
    if constexpr(R3D) {
        return FL(2.0) * REAL_PI;
    } else {
        return STD::sqrt(FL(2.0) * REAL_PI);
    }
}

// area of parallelogram formed by ray tube
template<typename MV> HOST_DEVICE inline real QScalar(const MV &q);
template<> HOST_DEVICE inline real QScalar(const vec2 &q) { return q.x; }
template<> HOST_DEVICE inline real QScalar(const mat2x2 &q)
{
    return glm::determinant(q);
}

template<bool R3D> HOST_DEVICE inline void ReceiverAngles(
    real &RcvrDeclAngle, real &RcvrAzimAngle, const VEC23<R3D> &rayt,
    const InfluenceRayInfo<R3D> &inflray)
{
    real angleXY = RadDeg * STD::atan2(rayt.y, rayt.x);
    if constexpr(R3D) {
        RcvrDeclAngle = RadDeg * STD::atan2(rayt.z, glm::length(XYCOMP(rayt)));
        RcvrAzimAngle = angleXY;
    } else {
        RcvrDeclAngle = angleXY;
        RcvrAzimAngle = inflray.init.SrcAzimAngle; // LP: Needed for Nx2D; for 2D, debug
                                                   // value
    }
}

////////////////////////////////////////////////////////////////////////////////
// Geom ray-centered helper functions
////////////////////////////////////////////////////////////////////////////////

template<bool R3D> HOST_DEVICE inline void RayNormalRayCen(
    const rayPt<R3D> &point, VEC23<R3D> &rayn1, VEC23<R3D> &rayn2)
{
    if constexpr(R3D) {
        RayNormal(point.t, point.phi, point.c, rayn1, rayn2);
    } else {
        // ray normal based on tangent with c(s) scaling
        rayn1 = vec2(point.t.y, -point.t.x) * point.c;
        rayn2 = vec2(DEBUG_LARGEVAL, DEBUG_LARGEVAL);
    }
}

HOST_DEVICE inline void Compute_RayCenCross(
    vec3 &xt, vec3 &xtxe1, vec3 &xtxe2, const vec3 &x, const vec3 &rayn1,
    const vec3 &rayn2, real zR, const InfluenceRayInfo<true> &inflray)
{
    // mbp: vector from receiver to each step of ray
    // LP: Not sure what this vector is, but it's not that.
    xt    = x - vec3(inflray.xs.x, inflray.xs.y, zR);
    xtxe1 = glm::cross(xt, rayn1);
    xtxe2 = glm::cross(xt, rayn2);
}

HOST_DEVICE inline void Compute_N_R_IR(
    real &n, real &r, int32_t &ir, const vec2 &x, const vec2 &normal, real zR,
    const Position *Pos)
{
    n = (zR - DEP(x)) / DEP(normal);
    r = x.x + n * normal.x;
    // following assumes uniform spacing in Pos->r
    ir = RToIR(r, Pos); // index of receiver
}

/**
 * *** Compute coordinates of intercept: nB, mB, rB ***
 */
HOST_DEVICE inline bool Compute_M_N_R_IR(
    real &m, real &n, real &r, int32_t &ir, const vec3 &e1xe2, const vec3 &xtxe1,
    const vec3 &xtxe2, const vec3 &xt, const vec2 &t_rcvr, const Position *Pos)
{
    real delta = -glm::dot(t_rcvr, XYCOMP(e1xe2));
    if(STD::abs(delta) < FL(1e-3)) return true; // ray normal || radial of rcvr line
    m  = glm::dot(t_rcvr, XYCOMP(xtxe1)) / delta;
    n  = -glm::dot(t_rcvr, XYCOMP(xtxe2)) / delta;
    r  = -glm::dot(xt, e1xe2) / delta;
    ir = RToIR(r, Pos); // index of nearest rcvr before normal
    return false;
}

////////////////////////////////////////////////////////////////////////////////
// Cerveny-only helper functions
////////////////////////////////////////////////////////////////////////////////

/**
 * Picks the optimum value for epsilon (LP: "beam constant" / "beam initial conditions")
 *
 * omega: angular frequency
 * c: sound speed
 * gradc: gradient
 * Dangle: angular spacing for ray fan (LP: either Dalpha or Dbeta)
 * rLoop: loop range (LP: in Beam struct)
 * EpsMultiplier: multiplier to manually adjust result (LP: in Beam struct)
 *
 * LP: For 3D, called twice, once for alpha and once for beta.
 * LP: Claims to support non-Cerveny cases, but never called or results never
 * used in non-Cerveny cases.
 */
template<bool O3D, bool R3D> HOST_DEVICE inline cpx PickEpsilon(
    real omega, real c, vec2 gradc, real angle, real Dangle,
    const BeamStructure<O3D> *Beam)
{
    // LP: BUG: Multiple codepaths do not set epsilonOpt, would lead to UB if
    // set up that way in env file. Instead, here they are set to debug values.
    real halfwidth        = R3D ? RL(0.0) : DEBUG_LARGEVAL;
    cpx epsilonOpt        = cpx(DEBUG_LARGEVAL, DEBUG_LARGEVAL);
    bool defaultHalfwidth = true, defaultEps = true, zeroEps = false;
    real cz;
    if(IsCervenyInfl(Beam)) {
        defaultHalfwidth = defaultEps = false;
        if(IsBeamWidthSpaceFilling(Beam)) {
            defaultHalfwidth = defaultEps = true;
        } else if(IsBeamWidthMinimum(Beam)) {
            halfwidth  = STD::sqrt(FL(2.0) * c * FL(1000.0) * Beam->rLoop / omega);
            defaultEps = true;
        } else if(IsBeamWidthWKB(Beam)) {
            if(R3D) {
                GlobalLog("Warning, WKB beams unimplemented in BELLHOP3D\n");
            } else {
                halfwidth  = REAL_MAX;
                cz         = gradc.y;
                epsilonOpt = (cz == FL(0.0))
                    ? RL(1e10)
                    : ((-STD::sin(angle) / STD::cos(SQ(angle))) * c * c / cz);
            }
        } else if(IsBeamWidthCerveny(Beam)) {
            GlobalLog(
                "Warning, Cerveny beam *width* (not type) %s\n",
                R3D ? "buggy in BELLHOP3D" : "not implemented in BELLHOP (2D)");
        } else {
            GlobalLog(
                "Invalid beam width option %c buggily ignored in PickEpsilon\n",
                Beam->Type[1]);
        }
    } else if(IsGeometricInfl(Beam)) {
        if constexpr(R3D) {
            zeroEps = true;
        } else {
            if(Beam->Type[0] == '^' || Beam->Type[0] == ' ') {
                GlobalLog("Warning, Beam->Type[0] = ^ or ' ' not properly handled in "
                          "BELLHOP (2D)\n");
                defaultHalfwidth = defaultEps = false;
            } else if(IsGaussianGeomInfl(Beam) && IsRayCenInfl(Beam)) {
                GlobalLog(BHC_PROGRAMNAME ": Geo Gaussian beams in ray-cent. coords. not "
                                          "implemented in BELLHOP (2D)\n");
                bail();
            }
        }
    } else if(IsSGBInfl(Beam)) {
        // LP: Supported here in 3D even though not supported in Influence 3D.
        NULLSTATEMENT;
    } else {
        GlobalLog("Invalid Beam->Type[0] %c ignored in PickEpsilon\n", Beam->Type[0]);
        defaultHalfwidth = defaultEps = false;
    }

    if(zeroEps) {
        epsilonOpt = FL(0.0);
    } else if(defaultEps) {
        if(defaultHalfwidth) {
            halfwidth = (Dangle == FL(0.0)) ? FL(0.0)
                                            : (FL(2.0) / ((omega / c) * Dangle));
        }
        epsilonOpt = J * FL(0.5) * omega * SQ(halfwidth);
    }

    /*
    // On first call write info to prt file
    static bool INIFlag = true;
    if(INIFlag){
        const char *tag = IsCervenyInfl(Beam) ? GetBeamWidthTag(Beam) :
    GetBeamTypeTag(Beam); PRTFile << "\n" << tag << "\n"; PRTFile << "halfwidth  = " <<
    halfwidth << "\n"; PRTFile << "epsilonOpt = " << epsilonOpt << "\n"; PRTFile <<
    "EpsMult    = " << Beam->epsMultiplier << "\n\n"; INIFlag = false;
    }
    */

    return Beam->epsMultiplier * epsilonOpt;
}

/**
 * Form gamma
 */
template<typename CFG, bool O3D> HOST_DEVICE inline cpx Compute_gamma(
    const rayPt<false> &point, const cpx &pB, const cpx &qB,
    const Origin<O3D, false> &org, const SSPStructure *ssp, SSPSegState &iSeg)
{
    vec2 rayt = point.c * point.t;     // unit tangent
    vec2 rayn = vec2(rayt.y, -rayt.x); // unit normal

    SSPOutputs<false> o;
    EvaluateSSP<CFG, O3D, false>(point.x, point.t, o, org, ssp, iSeg);

    real csq = SQ(o.ccpx.real());
    real cS  = glm::dot(o.gradc, rayt);
    real cN  = glm::dot(o.gradc, rayn);

    real Tr = rayt.x;
    real Tz = rayt.y;

    if(qB != RL(0.0)) {
        return FL(0.5)
            * (pB / qB * SQ(Tr) + FL(2.0) * cN / csq * Tz * Tr - cS / csq * SQ(Tz));
    } else {
        return FL(0.0);
    }
}

template<bool O3D> HOST_DEVICE inline void Compute_eps_pB_qB(
    cpx &eps, cpx &pB, cpx &qB, const rayPt<false> &point,
    const InfluenceRayInfo<false> &inflray, const BeamStructure<O3D> *Beam)
{
    if(IsBeamWidthCerveny(Beam)) {
        eps = J * STD::abs(point.q.x / point.q.y);
    } else {
        eps = inflray.epsilon1;
    }
    pB = point.p.x + eps * point.p.y;
    qB = point.q.x + eps * point.q.y;
}

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
    if(Ax <= x1) {
        return RL(1.0);
    } else if(Ax >= x2) {
        return RL(0.0);
    } else {
        real u = (Ax - x1) / (x2 - x1);
        return (RL(1.0) + RL(2.0) * u) * SQ(RL(1.0) - u);
    }
    // ret /= (FL(0.5) * (x1 + x2));
}

/**
 * Checks for a branch cut crossing and updates kmah accordingly
 */
template<bool O3D> HOST_DEVICE inline void BranchCut(
    const cpx &q1C, const cpx &q2C, int32_t &kmah, const BeamStructure<O3D> *Beam)
{
    real q1, q2;
    if(IsBeamWidthWKB(Beam)) {
        q1 = q1C.real();
        q2 = q2C.real();
    } else {
        if(q2C.real() >= FL(0.0)) return;
        q1 = q1C.imag();
        q2 = q2C.imag();
    }
    if((q1 < FL(0.0) && q2 >= FL(0.0)) || (q1 > FL(0.0) && q2 <= FL(0.0))) kmah = -kmah;
}

HOST_DEVICE inline vec2 FlipBeamForImage(
    const vec2 &x, int32_t image, const BdryType *Bdry)
{
    if(image == 1) { // True beam
        return x;
    } else if(image == 2) { // Surface-reflected beam
        return vec2(x.x, FL(2.0) * Bdry->Top.hs.Depth - x.y);
    } else if(image == 3) { // Bottom-reflected beam
        return vec2(x.x, FL(2.0) * Bdry->Bot.hs.Depth - x.y);
    } else {
        GlobalLog("Image index %d must be 1, 2, or 3\n", image);
        bail();
    }
}

HOST_DEVICE inline vec2 FlipNormalForImage(const vec2 &rayn, int32_t image)
{
    return vec2(rayn.x, rayn.y * ((image == 2) ? RL(-1.0) : RL(1.0)));
}

////////////////////////////////////////////////////////////////////////////////
// Caustics tracking system
////////////////////////////////////////////////////////////////////////////////

/**
 * LP: There are two versions of the phase shift condition used in the BELLHOP
 * code, with the equalities in opposite positions. qleq0 false is only used in
 * SGB.
 */
template<bool R3D> HOST_DEVICE inline bool IsAtCaustic(
    const InfluenceRayInfo<R3D> &inflray, real q, bool qleq0)
{
    if(qleq0) {
        return (q <= RL(0.0) && inflray.qOld > RL(0.0))
            || (q >= RL(0.0) && inflray.qOld < RL(0.0));
    } else {
        return (q < RL(0.0) && inflray.qOld >= RL(0.0))
            || (q > RL(0.0) && inflray.qOld <= RL(0.0));
    }
}

/**
 * phase shifts at caustics
 */
template<bool R3D> HOST_DEVICE inline void IncPhaseIfCaustic(
    InfluenceRayInfo<R3D> &inflray, real q, bool qleq0)
{
    if(IsAtCaustic<R3D>(inflray, q, qleq0)) inflray.phase += REAL_PI / FL(2.0);
}

/**
 * phase shifts at caustics
 */
template<bool R3D> HOST_DEVICE inline real FinalPhase(
    const rayPt<R3D> &point, const InfluenceRayInfo<R3D> &inflray, real q)
{
    real phaseInt = point.Phase + inflray.phase;
    if(IsAtCaustic<R3D>(inflray, q, true)) {
        // LP: All 2D influence functions discard point.Phase when this
        // condition is met. Probably a BUG as none of the 3D functions do this.
        phaseInt = (R3D ? phaseInt : inflray.phase) + REAL_PI / FL(2.0);
    }
    return phaseInt;
}

////////////////////////////////////////////////////////////////////////////////
// Storing results
////////////////////////////////////////////////////////////////////////////////

template<bool R3D> HOST_DEVICE inline void AddToField(
    cpxf *uAllSources, const cpxf &dfield, int32_t itheta, int32_t ir, int32_t iz,
    const InfluenceRayInfo<R3D> &inflray, const Position *Pos)
{
    size_t base = GetFieldAddr(
        inflray.init.isx, inflray.init.isy, inflray.init.isz, itheta, iz, ir, Pos);
    AtomicAddCpx(&uAllSources[base], dfield);
}

template<typename CFG, bool O3D, bool R3D> HOST_DEVICE inline void ApplyContribution(
    cpxf *uAllSources, real cnst, real w, real omega, cpx delay, real phaseInt,
    real RcvrDeclAngle, real RcvrAzimAngle, int32_t itheta, int32_t ir, int32_t iz,
    int32_t is, const InfluenceRayInfo<R3D> &inflray, const rayPt<R3D> &point1,
    const Position *Pos, const BeamStructure<O3D> *Beam, EigenInfo *eigen,
    const ArrInfo *arrinfo)
{
    if constexpr(O3D && !R3D) { itheta = inflray.init.ibeta; }
    if constexpr(CFG::run::IsEigenrays()) {
        // eigenrays
        RecordEigenHit(itheta, ir, iz, is, inflray.init, eigen);
    } else if constexpr(CFG::run::IsArrivals()) {
        // arrivals
        AddArr<R3D>(
            itheta, iz, ir, cnst * w, omega, phaseInt, delay, inflray.init, RcvrDeclAngle,
            RcvrAzimAngle, point1.NumTopBnc, point1.NumBotBnc, arrinfo, Pos);
    } else {
        cpxf dfield;
        if(IsCoherentRun(Beam)) {
            // coherent TL
            dfield = Cpx2Cpxf(cnst * w * STD::exp(-J * (omega * delay - phaseInt)));
            // omega * SQ(n) / (FL(2.0) * SQ(point1.c) * delay)))) // curvature correction
            // [LP: 2D only]
        } else {
            // incoherent/semicoherent TL
            real v = cnst * STD::exp((omega * delay).imag());
            v      = SQ(v) * w;
            if constexpr(CFG::infl::IsGaussianGeom()) {
                // Gaussian beam
                v *= GaussScaleFactor<R3D>();
            }
            dfield = cpxf((float)v, 0.0f);
        }
        AddToField<R3D>(uAllSources, dfield, itheta, ir, iz, inflray, Pos);
    }
}

////////////////////////////////////////////////////////////////////////////////
// Init_Influence
////////////////////////////////////////////////////////////////////////////////

template<bool O3D, bool R3D> inline void PreRun_Influence(const BeamStructure<O3D> *Beam)
{
    if(IsCervenyInfl(Beam)) {
        if constexpr(R3D) {
            // LP: The Influence3D (Cerveny) function is commented out; was
            // obviously implemented at some point.
            GlobalLog(
                IsCartesianInfl(Beam) ? "Run Type 'C' not supported at this time\n"
                                      : "Invalid Run Type\n");
            bail();
        } else if constexpr(O3D) {
            if(IsCartesianInfl(Beam)) {
#ifdef BHC_LIMIT_FEATURES
                GlobalLog("Nx2D Cerveny Cartesian is not supported by BELLHOP3D "
                          "but can be supported by " BHC_PROGRAMNAME " if you turn off "
                          "BHC_LIMIT_FEATURES\n");
                bail();
#else
                GlobalLog("Warning: Nx2D Cerveny Cartesian is not supported by "
                          "BELLHOP3D, but is supported by " BHC_PROGRAMNAME "\n");
#endif
            }
        }
    } else if(IsSGBInfl(Beam)) {
        if constexpr(R3D) {
            GlobalLog("Invalid Run Type\n");
            bail();
        }
    } else if(IsGeometricInfl(Beam)) {
        if constexpr(!R3D) {
            if(IsRayCenInfl(Beam) && IsGaussianGeomInfl(Beam)) {
#ifdef BHC_LIMIT_FEATURES
                GlobalLog("2D Gaussian RayCen is not supported by BELLHOP "
                          "but can be supported by " BHC_PROGRAMNAME " if you turn off "
                          "BHC_LIMIT_FEATURES\n");
                bail();
#else
                GlobalLog("Warning: 2D Gaussian RayCen is not supported by "
                          "BELLHOP, but is supported by " BHC_PROGRAMNAME "\n");
#endif
            }
        }
    } else {
        GlobalLog("Invalid Run Type\n");
        bail();
    }
}

template<typename CFG, bool O3D, bool R3D> HOST_DEVICE inline void Init_Influence(
    InfluenceRayInfo<R3D> &inflray, const rayPt<R3D> &point0, RayInitInfo &rinit,
    vec2 gradc, const Position *Pos, const Origin<O3D, R3D> &org, const SSPStructure *ssp,
    SSPSegState &iSeg, const AnglesStructure *Angles, const FreqInfo *freqinfo,
    const BeamStructure<O3D> *Beam)
{
    if constexpr(R3D) IGNORE_UNUSED(ssp);

    if constexpr(CFG::infl::IsCerveny() && !CFG::run::IsTL()) {
        GlobalLog("Cerveny influence does not support eigenrays or arrivals\n");
        bail();
    }

    inflray.init  = rinit;
    inflray.freq0 = freqinfo->freq0;
    inflray.omega = FL(2.0) * REAL_PI * inflray.freq0;
    inflray.c0    = point0.c;
    inflray.xs    = point0.x;
    // LP: The 5x version is changed to 50x on both codepaths before it is used.
    // inflray.RadMax = FL(5.0) * ccpx.real() / freqinfo->freq0; // 5 wavelength max
    // radius
    inflray.RadMax = FL(50.0) * point0.c / freqinfo->freq0; // 50 wavelength max radius
    inflray.Dalpha = Angles->alpha.d;
    inflray.Dbeta  = Angles->beta.d;

    const real BeamWindow = RL(4.0); // LP: Integer (!) in 2D
    inflray.BeamWindow    = CFG::infl::IsGaussianGeom() ? BeamWindow : RL(1.0);
    inflray.iBeamWindow2  = SQ(Beam->iBeamWindow);

    // LP: Did not have the abs for SGB, but it has been added.
    inflray.Ratio1 = STD::sqrt(STD::abs(STD::cos(rinit.alpha))); // point source
    if constexpr(R3D) {
        inflray.Ratio1 *= STD::sqrt(inflray.Dalpha * inflray.Dbeta) / point0.c;
    } else {
        if(IsLineSource(Beam)) inflray.Ratio1 = RL(1.0);
    }
    if constexpr(CFG::infl::IsGaussianGeom()) {
        inflray.Ratio1 /= GaussScaleFactor<R3D>();
    }

    if constexpr(CFG::infl::IsCerveny()) {
        inflray.epsilon1 = PickEpsilon<O3D, R3D>(
            inflray.omega, point0.c, gradc, rinit.alpha, Angles->alpha.d, Beam);
        if constexpr(R3D) {
            inflray.epsilon2 = PickEpsilon<O3D, R3D>(
                inflray.omega, point0.c, gradc, rinit.beta, Angles->beta.d, Beam);
        }
    }

    // LP: For all geometric types
    inflray.rcp_q0 = inflray.Dalpha / point0.c; // Reference for J = q0 / q = q * rcp_q0
    if constexpr(R3D) {
        inflray.rcp_qhat0 = STD::abs(STD::cos(inflray.init.alpha)) * inflray.Dbeta
            / point0.c;
    }

    // LP: For all geometric types
    inflray.phase = FL(0.0);

    if constexpr(CFG::infl::IsSGB()) {
        inflray.qOld = FL(1.0);
    } else {
        // LP: For Cartesian types
        inflray.qOld = QScalar(point0.q); // used to track KMAH index
    }

    // LP: For RayCen types
    if constexpr(CFG::infl::IsRayCen()) {
        if constexpr(CFG::infl::IsCerveny()) {
            // mbp: This logic means that the first step along the ray is skipped
            // which is a problem if deltas is very large, e.g. isospeed problems
            // I [mbp] fixed this in InfluenceGeoHatRayCen
            inflray.rayn1      = VEC23<R3D>(RL(0.0));
            DEP(inflray.rayn1) = RL(1.0);
            inflray.rayn2      = VEC23<R3D>(DEBUG_LARGEVAL);
            inflray.x          = VEC23<R3D>(RL(0.0));
            inflray.lastValid  = false;
        } else {
            RayNormalRayCen<R3D>(point0, inflray.rayn1, inflray.rayn2);
            inflray.x         = point0.x;
            inflray.lastValid = STD::abs(DEP(inflray.rayn1)) >= RL(1e-6); // LP: Ignored
                                                                          // in 3D
        }
    } else {
        // LP: For Cartesian types
        inflray.x = point0.x;
    }

    // LP: For Cerveny
    inflray.kmah = 1;

    if constexpr(CFG::infl::IsSGB()) {
        inflray.ir = 0;
    } else if constexpr(CFG::infl::IsGeometric() && CFG::infl::IsCartesian()) {
        // what if never satisfied?
        // what if there is a single receiver (ir = -1 possible)
        // LP: ir is always valid, even if it means not meeting the condition.
        real r;
        if constexpr(R3D) {
            // LP: Originally glm::distance(XYCOMP(point0.x), XYCOMP(inflray.xs))
            // but the ray always starts at the source
            r = RL(0.0);
        } else {
            r = point0.x.x;
        }
        inflray.ir = BinarySearchGT(Pos->Rr, Pos->NRr, 1, 0, r); // find index of first
                                                                 // receiver to the right
                                                                 // of rA
        if constexpr(!R3D) {
            if(point0.t.x < RL(0.0) && inflray.ir > 0)
                --inflray.ir; // if ray is left-traveling, get the first receiver to the
                              // left of rA
        }
    }

    if constexpr(CFG::infl::IsCerveny() && CFG::infl::IsCartesian()) {
        // LP: For Cerveny cart
        if constexpr(R3D) {
            GlobalLog("Run Type 'C' not supported at this time\n");
            bail();
        } else {
            // LP: Partially supported in Nx2D (O3D but not R3D)
            cpx eps0, pB0, qB0;
            Compute_eps_pB_qB<O3D>(eps0, pB0, qB0, point0, inflray, Beam);
            inflray.gamma = Compute_gamma<CFG, O3D>(point0, pB0, qB0, org, ssp, iSeg);
        }
    } else {
        // LP: not used
        inflray.gamma = RL(0.0);
    }
}

////////////////////////////////////////////////////////////////////////////////
// Cerveny step functions
////////////////////////////////////////////////////////////////////////////////

/**
 * Paraxial (Cerveny-style) beams in ray-centered coordinates
 */
template<typename CFG, bool O3D> HOST_DEVICE inline bool Step_InfluenceCervenyRayCen(
    const rayPt<false> &point0, const rayPt<false> &point1,
    InfluenceRayInfo<false> &inflray, int32_t is, cpxf *uAllSources, const BdryType *Bdry,
    const Position *Pos, const BeamStructure<O3D> *Beam)
{
    IGNORE_UNUSED(is);

    cpx eps0, eps1, pB0, pB1, qB0, qB1, gamma0, gamma1;
    // need to add logic related to NRz_per_range

    // During reflection imag(q) is constant and adjacent normals cannot bracket
    // a segment of the TL line, so no special treatment is necessary

    Compute_eps_pB_qB<O3D>(eps0, pB0, qB0, point0, inflray, Beam);
    Compute_eps_pB_qB<O3D>(eps1, pB1, qB1, point1, inflray, Beam);
    gamma0 = pB0 / qB0;
    gamma1 = pB1 / qB1;

    vec2 rayn1, dummy;
    RayNormalRayCen<false>(point1, rayn1, dummy);
    // If normal parallel to TL-line, skip to next step on ray
    // LP: Changed from (the FORTRAN equivalent of) REAL_MINPOS, see Fortran
    // version readme.
    if(STD::abs(DEP(rayn1)) < REAL_EPSILON) return true;

    // detect and skip duplicate points (happens at boundary reflection)
    if(IsDuplicatePoint(point0, point1, false)) {
        inflray.lastValid = true;
        inflray.x         = point1.x;
        inflray.rayn1     = rayn1;
        return true;
    }

    // compute KMAH index
    // Following is incorrect for 'Cerveny'-style beamwidth (narrow as possible)
    int32_t old_kmah = inflray.kmah;
    BranchCut<O3D>(qB0, qB1, inflray.kmah, Beam);

    for(int32_t iz = 0; iz < Pos->NRz_per_range; ++iz) {
        real zR = Pos->Rz[iz];

        for(int32_t image = 1; image <= Beam->Nimage; ++image) {
            real nA, rA, nB, rB;
            int32_t ir1, ir2; // LP: mbp switches from A/B naming to 1/2 here.
            Compute_N_R_IR(
                nB, rB, ir2, FlipBeamForImage(point1.x, image, Bdry),
                FlipNormalForImage(rayn1, image), zR, Pos);
            Compute_N_R_IR(
                nA, rA, ir1, FlipBeamForImage(inflray.x, image, Bdry),
                FlipNormalForImage(inflray.rayn1, image), zR, Pos);

            if(inflray.lastValid && ir1 < ir2) {
                for(int32_t ir = ir1 + 1; ir <= ir2; ++ir) {
                    real w, n, nSq, c;
                    cpx q, gamma, tau, contri;
                    w     = (Pos->Rr[ir] - rA) / (rB - rA);
                    q     = qB0 + w * (qB1 - qB0);
                    gamma = gamma0 + w * (gamma1 - gamma0);
                    n     = nA + w * (nB - nA);
                    nSq   = SQ(n);
                    if(gamma.imag() > 0) {
#ifndef BHC_USE_FLOATS
                        GlobalLog("Unbounded beam\n");
#endif
                        continue;
                    }

                    if(FL(-0.5) * inflray.omega * gamma.imag() * nSq
                       < inflray.iBeamWindow2) { // Within beam window?
                        c      = point0.c;
                        tau    = point0.tau + w * (point1.tau - point0.tau);
                        contri = inflray.Ratio1 * point1.Amp
                            * STD::sqrt(c * STD::abs(eps1) / q)
                            * STD::exp(
                                     -J
                                     * (inflray.omega * (tau + FL(0.5) * gamma * nSq)
                                        - point1.Phase));

                        cpx P_n = -J * inflray.omega * gamma * n * contri;
                        cpx P_s = -J * inflray.omega / c * contri;
                        switch(Beam->Component) {
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
                        BranchCut<O3D>(qB0, q, kmah, Beam); // Get correct branch of
                                                            // STD::sqrt

                        if(kmah < 0) contri = -contri;
                        if(image == 2) contri = -contri;

                        // LP: Non-TL runs not handled correctly.
                        if constexpr(CFG::run::IsTL()) {
                            if(!IsCoherentRun(Beam)) {
                                contri = contri * STD::conj(contri);
                            }
                        }
                        contri *= Hermite(n, inflray.RadMax, FL(2.0) * inflray.RadMax);

                        AddToField<false>(
                            uAllSources, Cpx2Cpxf(contri), O3D ? inflray.init.ibeta : 0,
                            ir, iz, inflray, Pos);
                    }
                }
            }
        }
    }

    inflray.lastValid = true;
    inflray.x         = point1.x;
    inflray.rayn1     = rayn1;
    return true;
}

/**
 * Paraxial (Cerveny-style) beams in Cartesian coordinates
 */
template<typename CFG, bool O3D> HOST_DEVICE inline bool Step_InfluenceCervenyCart(
    const rayPt<false> &point0, const rayPt<false> &point1,
    InfluenceRayInfo<false> &inflray, int32_t is, cpxf *uAllSources, const BdryType *Bdry,
    const Origin<O3D, false> &org, const SSPStructure *ssp, SSPSegState &iSeg,
    const Position *Pos, const BeamStructure<O3D> *Beam)
{
    cpx eps0, eps1, pB0, pB1, qB0, qB1, gamma0, gamma1;
    real zR;
    // need to add logic related to NRz_per_range

    // During reflection imag(q) is constant and adjacent normals cannot bracket
    // a segment of the TL line, so no special treatment is necessary

    Compute_eps_pB_qB<O3D>(eps0, pB0, qB0, point0, inflray, Beam);
    Compute_eps_pB_qB<O3D>(eps1, pB1, qB1, point1, inflray, Beam);

    // Form gamma and KMAH index
    // Treatment of KMAH index is incorrect for 'Cerveny' style beam width BeamType
    gamma0        = inflray.gamma;
    gamma1        = Compute_gamma<CFG, O3D>(point1, pB1, qB1, org, ssp, iSeg);
    inflray.gamma = gamma1;

    int32_t old_kmah = inflray.kmah;
    BranchCut<O3D>(qB0, qB1, inflray.kmah, Beam);

    if(is == 0) return true; // LP: Skips the first valid pair.
    // LP: Assumes rays may never travel left.
    if(point1.x.x > Pos->Rr[Pos->NRr - 1]) {
        return false; // LP: Terminates ray.
    }
    real rA = point0.x.x;
    real rB = point1.x.x;
    if(IsDuplicatePoint(point0, point1, false))
        return true; // don't process duplicate points

    // Compute upper index on rcvr line
    // Assumes r is a vector of equally spaced points
    int32_t irA = RToIR(rA, Pos);
    int32_t irB = RToIR(rB, Pos);

    if(irA >= irB) return true;

    for(int32_t ir = irA + 1; ir <= irB; ++ir) {
        real w, c;
        vec2 x, rayt;
        cpx q, tau, gamma, cnst;
        w     = (Pos->Rr[ir] - rA) / (rB - rA);
        x     = point0.x + w * (point1.x - point0.x);
        rayt  = point0.t + w * (point1.t - point0.t);
        c     = point0.c + w * (point1.c - point0.c);
        q     = qB0 + w * (qB1 - qB0);
        tau   = point0.tau + w * (point1.tau - point0.tau);
        gamma = gamma0 + w * (gamma1 - gamma0);

        if(gamma.imag() > FL(0.0)) {
#ifndef BHC_USE_FLOATS
            GlobalLog("Unbounded beam\n");
#endif
            continue;
        }

        cnst = inflray.Ratio1 * STD::sqrt(c * STD::abs(eps0) / q);

        // Get correct branch of STD::sqrt
        int32_t kmah = old_kmah;
        BranchCut<O3D>(qB0, q, kmah, Beam);
        if(kmah < 0) cnst = -cnst;

        for(int32_t iz = 0; iz < Pos->NRz_per_range; ++iz) {
            zR = Pos->Rz[iz];

            cpx contri = FL(0.0);
            for(int32_t image = 1; image <= Beam->Nimage; ++image) {
                real deltaz, Polarity;
                if(image == 1) { // True beam
                    deltaz   = zR - x.y;
                    Polarity = RL(1.0);
                } else if(image == 2) { // Surface reflected beam
                    deltaz   = -zR + FL(2.0) * Bdry->Top.hs.Depth - x.y;
                    Polarity = RL(-1.0);
                } else if(image == 3) { // Bottom  reflected beam
                    deltaz   = -zR + FL(2.0) * Bdry->Bot.hs.Depth - x.y;
                    Polarity = RL(1.0); // assumes rigid bottom
                } else {
                    GlobalLog("Invalid Beam->Nimage %d\n", Beam->Nimage);
                    bail();
                }
                if(inflray.omega * gamma.imag() * SQ(deltaz) < inflray.iBeamWindow2) {
                    contri += Polarity * point1.Amp
                        * Hermite(deltaz, inflray.RadMax, FL(2.0) * inflray.RadMax)
                        * STD::exp(
                                  -J
                                  * (inflray.omega
                                         * (tau + rayt.y * deltaz + gamma * SQ(deltaz))
                                     - point1.Phase));
                }
            }

            // LP: Non-TL not handled correctly.
            if constexpr(CFG::run::IsTL()) {
                contri = cnst * contri;
                if(!IsCoherentRun(Beam)) { contri = contri * STD::conj(contri); }
            }

            AddToField<false>(
                uAllSources, Cpx2Cpxf(contri), O3D ? inflray.init.ibeta : 0, ir, iz,
                inflray, Pos);
        }
    }

    return true;
}

////////////////////////////////////////////////////////////////////////////////
// Geom step functions
////////////////////////////////////////////////////////////////////////////////

/**
 * LP: Core influence function for all geometrical types: 2D/3D, Cartesian /
 * ray-centered, hat / Gaussian
 */
template<typename CFG, bool O3D, bool R3D> HOST_DEVICE inline void InfluenceGeoCore(
    real s, real n1, real n2, const V2M2<R3D> &dq, const cpx &dtau, int32_t itheta,
    int32_t ir, int32_t iz, int32_t is, const rayPt<R3D> &point0,
    const rayPt<R3D> &point1, real RcvrDeclAngle, real RcvrAzimAngle,
    const InfluenceRayInfo<R3D> &inflray, cpxf *uAllSources, const Position *Pos,
    const BeamStructure<O3D> *Beam, EigenInfo *eigen, const ArrInfo *arrinfo)
{
    static_assert(
        CFG::infl::IsGeometric(), "InfluenceGeoCore templated with non-geometric type!");

    V2M2<R3D> qInterp = point0.q + s * dq; // interpolated amplitude
    if constexpr(R3D) {
        qInterp[0][0] *= inflray.rcp_q0;
        qInterp[1][0] *= inflray.rcp_q0;
        qInterp[0][1] *= inflray.rcp_qhat0;
        qInterp[1][1] *= inflray.rcp_qhat0;
    }
    real qFinal = QScalar(qInterp);

    real n1prime, n2prime;  // LP: equal to n1 in 2D, rotated & normalized in 3D (were a,
                            // b)
    real sigma, sigma_orig; // beam radius
    real beamCoordDist;
    if constexpr(R3D) {
        if constexpr(CFG::infl::IsCartesian()) {
            real l1 = glm::length(glm::row(qInterp, 0));
            real l2 = glm::length(glm::row(qInterp, 1));
            if(l1 == FL(0.0) || l2 == FL(0.0)) {
#ifdef INFL_DEBUGGING_IR
                if(itheta == INFL_DEBUGGING_ITHETA && ir == INFL_DEBUGGING_IR
                   && iz == INFL_DEBUGGING_IZ) {
                    GlobalLog(
                        "is itheta iz ir %3d %3d %3d %3d: Skipping b/c l1/l2\n", is,
                        itheta, iz, ir);
                }
#endif
                return;
            }
        }

        if(qFinal == FL(0.0)) {
// receiver is outside the beam
#ifdef INFL_DEBUGGING_IR
            if(itheta == INFL_DEBUGGING_ITHETA && ir == INFL_DEBUGGING_IR
               && iz == INFL_DEBUGGING_IZ) {
                GlobalLog(
                    "is itheta iz ir %3d %3d %3d %3d: Skipping b/c qFinal\n", is, itheta,
                    iz, ir);
            }
#endif
            return;
        }

        n1prime = STD::abs((-qInterp[0][1] * n2 + qInterp[1][1] * n1) / qFinal);
        n2prime = STD::abs((qInterp[0][0] * n2 - qInterp[1][0] * n1) / qFinal);
#ifdef INFL_DEBUGGING_IR
        if(itheta == INFL_DEBUGGING_ITHETA && ir == INFL_DEBUGGING_IR
           && iz == INFL_DEBUGGING_IZ) {
            GlobalLog("is itheta iz ir %3d %3d %3d %3d:\n", is, itheta, iz, ir);
            PrintMatrix(qInterp, "qInterp");
            GlobalLog("qFinal n1 n2 %13.8g %13.8g %13.8g\n", qFinal, n1, n2);
            GlobalLog("n1prime %13.8f n2prime %13.8f\n", n1prime, n2prime);
        }
#endif

        if constexpr(CFG::infl::IsGaussianGeom()) {
            beamCoordDist = n1prime + n2prime;
        } else {
            beamCoordDist = bhc::max(n1prime, n2prime);
        }
        sigma = FL(1.0);
    } else {
        sigma = sigma_orig = STD::abs(qFinal * inflray.rcp_q0); // LP: called RadiusMax or
                                                                // l in non-Gaussian
        AdjustSigma<CFG, O3D, R3D>(sigma, point0, point1, inflray);

        beamCoordDist = n1prime = n1;
        IGNORE_UNUSED(n2);
        IGNORE_UNUSED(n2prime);
    }
#ifdef INFL_DEBUGGING_IR
    if(itheta == INFL_DEBUGGING_ITHETA && ir == INFL_DEBUGGING_IR
       && iz == INFL_DEBUGGING_IZ) {
        GlobalLog(
            "is itheta iz ir %3d %3d %3d %3d: (a+b) %g BeamWindow %g\n", is, itheta, iz,
            ir, beamCoordDist, inflray.BeamWindow * sigma);
    }
#endif
    if((beamCoordDist > inflray.BeamWindow * sigma ||
        // LP: 2D versions use >= (or rather < for write)
        (!R3D && beamCoordDist == inflray.BeamWindow * sigma))
       // LP: The "outside of beam window" condition is commented out in 3D Gaussian
       // raycen.
       && !(R3D && CFG::infl::IsRayCen() && CFG::infl::IsGaussianGeom())) {
        // receiver is outside the beam
        return;
    }

    cpx delay    = point0.tau + s * dtau; // interpolated delay
    real cfactor = point1.c;
    if constexpr(!R3D) cfactor = STD::sqrt(cfactor);
    real cnst = inflray.Ratio1 * cfactor * point1.Amp / STD::sqrt(STD::abs(qFinal));
    real w;
    if constexpr(R3D) {
        if constexpr(CFG::infl::IsGaussianGeom()) {
            w = STD::exp(FL(-0.5) * (SQ(n1prime) + SQ(n2prime)));
        } else {
            w = (FL(1.0) - n1prime) * (FL(1.0) - n2prime);
        }
    } else {
        if constexpr(CFG::infl::IsGaussianGeom()) {
            w = STD::exp(FL(-0.5) * SQ(n1prime / sigma))
                * (sigma_orig / sigma); // Gaussian
                                        // decay
        } else {
            w = (sigma - n1prime) / sigma; // hat function: 1 on center, 0 on edge
        }
    }
    real phaseInt = FinalPhase<R3D>(
        (!R3D && CFG::infl::IsGaussianGeom() ? point1 : point0), inflray, qFinal);

#ifdef INFL_DEBUGGING_IR
    if(itheta == INFL_DEBUGGING_ITHETA && ir == INFL_DEBUGGING_IR
       && iz == INFL_DEBUGGING_IZ) {
        GlobalLog(
            "is itheta iz ir %3d %3d %3d %3d: cnst w delay phaseInt %g %g (%g,%g) %g\n",
            is, itheta, iz, ir, cnst, w, delay.real(), delay.imag(), phaseInt);
    }
#endif

    ApplyContribution<CFG, O3D, R3D>(
        uAllSources, cnst, w, inflray.omega, delay, phaseInt, RcvrDeclAngle,
        RcvrAzimAngle, itheta, ir, iz, is, inflray, point1, Pos, Beam, eigen, arrinfo);
}

/**
 * Geometrically-spreading beams with a hat- or Gaussian-shaped beam
 * in ray-centered coordinates
 */
template<typename CFG, bool O3D, bool R3D> HOST_DEVICE inline bool
Step_InfluenceGeoRayCen(
    const rayPt<R3D> &point0, const rayPt<R3D> &point1, InfluenceRayInfo<R3D> &inflray,
    int32_t is, cpxf *uAllSources, const Position *Pos, const BeamStructure<O3D> *Beam,
    EigenInfo *eigen, const ArrInfo *arrinfo)
{
    real phaseq = QScalar(point0.q);
    IncPhaseIfCaustic<R3D>(inflray, phaseq, true);
    inflray.qOld = phaseq;

    V2M2<R3D> dq = point1.q - point0.q;
    cpx dtau     = point1.tau - point0.tau;

    vec3 e1xe2A, e1xe2B;
    if constexpr(R3D) {
        e1xe2A = point0.c * point0.t;
        e1xe2B = point1.c * point1.t;
    } else {
        IGNORE_UNUSED(e1xe2A);
        IGNORE_UNUSED(e1xe2B);
    }

    VEC23<R3D> rayn1, rayn2; // LP: rayn1 was rn, zn in 2D
    RayNormalRayCen<R3D>(point1, rayn1, rayn2);
    if constexpr(!R3D) {
        if(STD::abs(DEP(rayn1)) < RL(1e-10)) return true;
    }

    // LP: In 3D hat, non-constant threshold (same as all others) commented out
    // in favor of constant one.
    if(IsDuplicatePoint(point0, point1, R3D && CFG::infl::IsHatGeom())) {
        inflray.lastValid = true;
        inflray.x         = point1.x;
        inflray.rayn1     = rayn1;
        inflray.rayn2     = rayn2;
#ifdef INFL_DEBUGGING_IR
        GlobalLog("Skipping b/c dupl\n");
#endif
        return true;
    }

    real RcvrDeclAngle, RcvrAzimAngle;
    ReceiverAngles<R3D>(
        RcvrDeclAngle, RcvrAzimAngle, R3D ? (point1.x - point0.x) : point1.t, inflray);

    // During reflection imag(q) is constant and adjacent normals cannot bracket
    // a segment of the TL line, so no special treatment is necessary

    for(int32_t iz = 0; iz < Pos->NRz_per_range; ++iz) {
        real zR = Pos->Rz[iz];

        vec3 xtA, xtB, xtxe1A, xtxe1B, xtxe2A, xtxe2B;
        if constexpr(R3D) {
            Compute_RayCenCross(
                xtA, xtxe1A, xtxe2A, point0.x, inflray.rayn1, inflray.rayn2, zR, inflray);
            Compute_RayCenCross(xtB, xtxe1B, xtxe2B, point1.x, rayn1, rayn2, zR, inflray);
        } else {
            IGNORE_UNUSED(xtA);
            IGNORE_UNUSED(xtxe1A);
            IGNORE_UNUSED(xtxe2A);
            IGNORE_UNUSED(xtB);
            IGNORE_UNUSED(xtxe1B);
            IGNORE_UNUSED(xtxe2B);
        }

        int32_t itheta = -1;
        do {
            ++itheta;
            if constexpr(R3D) {
                if(itheta >= Pos->Ntheta) break;
            }

            real mA, mB, nA, nB, rA, rB;
            int32_t irA, irB;
            if constexpr(R3D) {
                if(Compute_M_N_R_IR(
                       mA, nA, rA, irA, e1xe2A, xtxe1A, xtxe2A, xtA, Pos->t_rcvr[itheta],
                       Pos)) {
#ifdef INFL_DEBUGGING_IR
                    if(itheta == INFL_DEBUGGING_ITHETA && iz == INFL_DEBUGGING_IZ) {
                        GlobalLog("Skipping b/c deltaA\n");
                    }
#endif
                    continue;
                }
                if(Compute_M_N_R_IR(
                       mB, nB, rB, irB, e1xe2B, xtxe1B, xtxe2B, xtB, Pos->t_rcvr[itheta],
                       Pos)) {
#ifdef INFL_DEBUGGING_IR
                    if(itheta == INFL_DEBUGGING_ITHETA && iz == INFL_DEBUGGING_IZ) {
                        GlobalLog("Skipping b/c deltaB\n");
                    }
#endif
                    continue;
                }
                if constexpr(CFG::infl::IsGaussianGeom()) {
                    // Possible contribution if max possible beamwidth > min possible
                    // distance to receiver
                    real MaxRadius_m = inflray.BeamWindow
                        * bhc::max(STD::abs(point0.q[1][1] * inflray.rcp_qhat0),
                                   STD::abs(point1.q[1][1] * inflray.rcp_qhat0));
                    real MaxRadius_n = inflray.BeamWindow
                        * bhc::max(STD::abs(point0.q[0][0] * inflray.rcp_q0),
                                   STD::abs(point1.q[0][0] * inflray.rcp_q0));
                    if(MaxRadius_m <= bhc::min(STD::abs(mA), STD::abs(mB))
                       && mA * mB >= RL(0.0)) {
#ifdef INFL_DEBUGGING_IR
                        if(itheta == INFL_DEBUGGING_ITHETA && iz == INFL_DEBUGGING_IZ) {
                            GlobalLog(
                                "Skipping b/c MaxRadius_m %g mA %g mB %g\n", MaxRadius_m,
                                mA, mB);
                        }
#endif
                        continue;
                    }
                    if(MaxRadius_n <= bhc::min(STD::abs(nA), STD::abs(nB))
                       && nA * nB >= RL(0.0)) {
#ifdef INFL_DEBUGGING_IR
                        if(itheta == INFL_DEBUGGING_ITHETA && iz == INFL_DEBUGGING_IZ) {
                            GlobalLog(
                                "Skipping b/c MaxRadius_n %g nA %g nB %g\n", MaxRadius_n,
                                nA, nB);
                        }
#endif
                        continue;
                    }
                }
            } else {
                if(!inflray.lastValid) {
                    nA  = RL(1e10);
                    rA  = RL(1e10);
                    irA = 0;
                } else {
                    Compute_N_R_IR(nA, rA, irA, inflray.x, inflray.rayn1, zR, Pos);
                }
                Compute_N_R_IR(nB, rB, irB, point1.x, rayn1, zR, Pos);
            }

            if(irA == irB) {
#ifdef INFL_DEBUGGING_IR
                if(itheta == INFL_DEBUGGING_ITHETA && iz == INFL_DEBUGGING_IZ) {
                    GlobalLog("Skipping b/c irA == irB\n");
                }
#endif
                continue;
            }

            // *** Compute contributions to bracketed receivers ***

            for(int32_t ir = bhc::min(irA, irB) + 1; ir <= bhc::max(irA, irB); ++ir) {
                real s  = (Pos->Rr[ir] - rA) / (rB - rA); // LP: called w in 2D
                real n1 = STD::abs(nA + s * (nB - nA));   // normal distance to ray
                real n2 = DEBUG_LARGEVAL; // LP: n1, n2 were called n, m in 3D
                if constexpr(R3D) {
                    n2 = STD::abs(mA + s * (mB - mA)); // normal distance to ray
                }
                InfluenceGeoCore<CFG, O3D, R3D>(
                    s, n1, n2, dq, dtau, itheta, ir, iz, is, point0, point1,
                    RcvrDeclAngle, RcvrAzimAngle, inflray, uAllSources, Pos, Beam, eigen,
                    arrinfo);
            }
        } while(R3D);
    }

    inflray.lastValid = true;
    inflray.x         = point1.x;
    inflray.rayn1     = rayn1;
    inflray.rayn2     = rayn2;
    return true;
}

/**
 * Geometric, hat-shaped or Gaussian beams in Cartesian coordintes
 *
 * uAllSources: complex pressure field
 */
template<typename CFG, bool O3D, bool R3D> HOST_DEVICE inline bool Step_InfluenceGeoCart(
    const rayPt<R3D> &point0, const rayPt<R3D> &point1, InfluenceRayInfo<R3D> &inflray,
    int32_t is, cpxf *uAllSources, const Position *Pos, const BeamStructure<O3D> *Beam,
    EigenInfo *eigen, const ArrInfo *arrinfo)
{
    // LP: Replaced ScaleBeam in 3D with applying the same scale factors below.
    // This avoids modifying the ray and makes the codepaths more similar.

    real rA, rB;
    if constexpr(R3D) {
        // LP: 3D updates rA even in the early return conditions below, so
        // it doesn't use inflray.x at all.
        rA = glm::length(XYCOMP(point0.x) - XYCOMP(inflray.xs));
        rB = glm::length(XYCOMP(point1.x) - XYCOMP(inflray.xs));
        if(IsDuplicatePoint<R3D>(rB, rA, false)) return true;
        // LP: This is silly logic, see comments in influence3D.f90
        if(is == 0) { inflray.ir = (rB > rA) ? 0 : Pos->NRr - 1; }
    } else {
        // LP: This is different from point0.x.x due to early return for duplicate points.
        rA = inflray.x.x;
        rB = point1.x.x;
    }
    VEC23<R3D> x_ray = point0.x;

    // compute normalized tangent (compute it because we need to measure the step length)
    VEC23<R3D> rayt = point1.x - point0.x;
    real rlen       = glm::length(rayt);
    // if duplicate point in ray, skip to next step along the ray
    // LP: 2D: and don't update rA (inflray.x) for next time
    if(IsSmallValue<R3D>(rlen, point1.x.x, false)) return true;
    rayt /= rlen;

    // LP: Ray normals
    VEC23<R3D> rayn1, rayn2; // LP: e1, e2 in 3D; rayn, (none) in 2D
    if constexpr(R3D) {
        RayNormal_unit(rayt, point1.phi, rayn1, rayn2);
    } else {
        rayn1 = vec2(-rayt.y, rayt.x); // unit normal to ray
        IGNORE_UNUSED(rayn2);
    }
    real RcvrDeclAngle, RcvrAzimAngle;
    ReceiverAngles<R3D>(RcvrDeclAngle, RcvrAzimAngle, rayt, inflray);

    // LP: Quantities to be interpolated between steps
    V2M2<R3D> dq = point1.q - point0.q;     // LP: dqds in 2D
    cpx dtau     = point1.tau - point0.tau; // LP: dtauds in 2D

    // phase shifts at caustics
    real phaseq = QScalar(point0.q);
    IncPhaseIfCaustic<R3D>(inflray, phaseq, true);
    inflray.qOld = phaseq;

    real L_diag;     // LP: 3D
    real zmin, zmax; // LP: 2D here, 3D later
    // beam window: kills beams outside exp(RL(-0.5) * SQ(ibwin))
    if constexpr(R3D) {
        // beamwidths / LP: Variable values don't carry over to per-receiver beamwidth
        // below
        real l1 = bhc::max(
            glm::length(glm::row(point0.q, 0) * inflray.rcp_q0),
            glm::length(glm::row(point1.q, 0) * inflray.rcp_q0));
        real l2 = bhc::max(
            glm::length(glm::row(point0.q, 1) * inflray.rcp_qhat0),
            glm::length(glm::row(point1.q, 1) * inflray.rcp_qhat0));
        // worst case is when rectangle is rotated to catch the hypotenuse
        L_diag = STD::sqrt(SQ(l1) + SQ(l2));
    } else {
        real sigma, RadiusMax;
        sigma = bhc::max(STD::abs(point0.q.x), STD::abs(point1.q.x)) * inflray.rcp_q0
            / STD::abs(rayt.x); // beam radius projected onto vertical line
        AdjustSigma<CFG, O3D, R3D>(sigma, point0, point1, inflray);
        RadiusMax = inflray.BeamWindow * sigma; // LP: 1 * sigma for non-Gaussian
        // depth limits of beam
        // LP: For rays shot at exactly 60 degrees, they will hit this edge case.
        // This is a sharp edge--the handling on each side of this edge may be
        // significantly different. So, moved the edge away from the round number.
        if(STD::abs(rayt.x) > FL(0.50001)) { // shallow angle ray
            zmin = bhc::min(point0.x.y, point1.x.y) - RadiusMax;
            zmax = bhc::max(point0.x.y, point1.x.y) + RadiusMax;
        } else { // steep angle ray
            zmin = -REAL_MAX;
            zmax = REAL_MAX;
        }
        IGNORE_UNUSED(L_diag);
    }

    // compute beam influence for this segment of the ray
    while(true) {
        // is Rr[ir] contained in [rA, rB)? Then compute beam influence
        // LP: Because of the new setup and always incrementing regardless of
        // which direction the ray goes, we only have to check this side.
        if(Pos->Rr[inflray.ir] >= bhc::min(rA, rB)
           && Pos->Rr[inflray.ir] < bhc::max(rA, rB)) {
            int32_t itheta = -1;
            do {
                VEC23<R3D> x_rcvr;

                ++itheta;
                if constexpr(R3D) {
                    // LP: Loop logic
                    if(itheta >= Pos->Ntheta) break;

                    vec2 t_rcvr = Pos->t_rcvr[itheta];
                    SETXY(
                        x_rcvr, XYCOMP(inflray.xs) + (real)Pos->Rr[inflray.ir] * t_rcvr);
                    // normal distance from rcvr to ray segment
                    real m_prime = STD::abs(
                        glm::dot(XYCOMP(x_rcvr) - XYCOMP(x_ray), vec2(-rayt.y, rayt.x)));
                    // LP: Commented out in Gaussian
                    if(CFG::infl::IsHatGeom() && m_prime > inflray.BeamWindow * L_diag) {
                        // GlobalLog("Skip theta b/c m_prime %g > L_diag %g\n", m_prime,
                        // L_diag);
                        continue;
                    }

                    // The set of possible receivers is a ring
                    // However, extrapolating the beam backwards produces contributions
                    // with s negative and large We do not want to accept these
                    // contributions--- they have the proper range but are 180 degrees
                    // away from this segment of the ray
                    // LP: This s value is not reused below
                    // vector to rcvr dotted into vector to ray point
                    real s = glm::dot(
                        XYCOMP(x_rcvr) - XYCOMP(inflray.xs),
                        XYCOMP(x_ray) - XYCOMP(inflray.xs));
                    if(s < RL(0.0)) {
                        // GlobalLog("Skip theta b/c s\n");
                        continue;
                    }

                    // calculate z-limits for the beam (could be pre-cacluated for each
                    // itheta)
                    // normal to the vertical receiver plane
                    vec2 e_theta = vec2(-t_rcvr.y, t_rcvr.x);
                    // normal to the ray in the vertical receiver plane
                    real n_ray_z = rayt.x * e_theta.y - rayt.y * e_theta.x;

                    if(STD::abs(n_ray_z) < RL(1e-9)) {
                        // GlobalLog("Skip theta b/c L_z divide by zero %g\n", n_ray_z);
                        continue; // avoid divide by zero
                    }
                    real L_z = inflray.BeamWindow * L_diag / STD::abs(n_ray_z);
                    // min depth of ray segment
                    zmin = bhc::min(point0.x.z, point1.x.z) - L_z;
                    // max depth of ray segment
                    zmax = bhc::max(point0.x.z, point1.x.z) + L_z;

                    // if(inflray.ir == 7 && itheta == 59){
                    //     GlobalLog("step ir itheta %d %d %d\n", is, inflray.ir, itheta);
                    // }
                } else {
                    IGNORE_UNUSED(itheta);
                    x_rcvr.x = Pos->Rr[inflray.ir];
                }

                for(int32_t iz = 0; iz < Pos->NRz_per_range; ++iz) {
                    int32_t tempiz = iz;
                    if constexpr(!R3D) {
                        if(IsIrregularGrid(Beam)) tempiz = inflray.ir;
                    } // else rectilinear grid
                    DEP(x_rcvr) = Pos->Rz[tempiz];
                    if(DEP(x_rcvr) < zmin || DEP(x_rcvr) > zmax) continue;

                    VEC23<R3D> x_rcvr_ray = x_rcvr - x_ray;

                    // linear interpolation of q's
                    real s = glm::dot(x_rcvr_ray, rayt) / rlen; // proportional distance
                                                                // along ray
                    real n1 = STD::abs(glm::dot(x_rcvr_ray, rayn1)); // normal distance to
                                                                     // ray
                    real n2 = DEBUG_LARGEVAL;
                    if constexpr(R3D) {
                        n2 = STD::abs(glm::dot(x_rcvr_ray, rayn2)); // normal distance to
                                                                    // ray
                    }
                    InfluenceGeoCore<CFG, O3D, R3D>(
                        s, n1, n2, dq, dtau, itheta, inflray.ir, iz, is, point0, point1,
                        RcvrDeclAngle, RcvrAzimAngle, inflray, uAllSources, Pos, Beam,
                        eigen, arrinfo);
                }
            } while(R3D);
        }

        // receiver not bracketed; bump receiver index, inflray.ir, towards rB
        int32_t irTT;
        if(rB > Pos->Rr[inflray.ir]) {
            if(inflray.ir >= Pos->NRr - 1) break; // go to next step on ray
            irTT = inflray.ir + 1;                // bump right
            if(Pos->Rr[irTT] >= rB) break;        // go to next step on ray
        } else {
            if(inflray.ir <= 0) break;     // go to next step on ray
            irTT = inflray.ir - 1;         // bump left
            if(Pos->Rr[irTT] <= rB) break; // go to next step on ray
        }
        inflray.ir = irTT;
    }

    inflray.x = point1.x;
    return true;
}

////////////////////////////////////////////////////////////////////////////////
// Other step function
////////////////////////////////////////////////////////////////////////////////

/**
 * Bucker's Simple Gaussian Beams in Cartesian coordinates
 */
template<typename CFG, bool O3D> HOST_DEVICE inline bool Step_InfluenceSGB(
    const rayPt<false> &point0, const rayPt<false> &point1,
    InfluenceRayInfo<false> &inflray, int32_t is, cpxf *uAllSources, const Position *Pos,
    const BeamStructure<O3D> *Beam, EigenInfo *eigen, const ArrInfo *arrinfo)
{
    real w;
    vec2 x, rayt;
    cpx tau;
    real RcvrDeclAngle, RcvrAzimAngle;
    ReceiverAngles<false>(RcvrDeclAngle, RcvrAzimAngle, point1.t, inflray);

    const real beta = FL(0.98); // Beam Factor
    real a          = FL(-4.0) * STD::log(beta) / SQ(inflray.Dalpha);
    real cn         = inflray.Dalpha * STD::sqrt(a / REAL_PI);
    real rA         = point0.x.x;

    real rB = point1.x.x;

    real q = point0.q.x;
    IncPhaseIfCaustic<false>(inflray, q, false);
    inflray.qOld = q;

    // Loop over bracketed receiver ranges
    // LP: BUG: See README.md.
    while(STD::abs(rB - rA) > RL(1.0e3) * spacing(rA) && rB > Pos->Rr[inflray.ir]) {
        w    = (Pos->Rr[inflray.ir] - rA) / (rB - rA);
        x    = point0.x + w * (point1.x - point0.x);
        rayt = point0.t + w * (point1.t - point0.t);
        q    = point0.q.x + w * (point1.q.x - point0.q.x);
        tau  = point0.tau + w * (point1.tau - point0.tau);

        // following is incorrect because ray doesn't always use a step of deltas
        // LP: The while ignores extremely small steps, but those small steps
        // still increment is, so the later ray segments still treat it as if
        // all steps leading up to them were of size deltas.
        real sint = ((real)(is + 1) + w) * Beam->deltas;

        IncPhaseIfCaustic<false>(inflray, q, false);

        // GlobalLog("is ir %d %d\n", is, inflray.ir);

        for(int32_t iz = 0; iz < Pos->NRz_per_range; ++iz) {
            real deltaz = Pos->Rz[iz] - x.y; // ray to rcvr distance
            // LP: Reinstated this condition for eigenrays and arrivals, as
            // without it every ray would be an eigenray / arrival.
            real Adeltaz = STD::abs(deltaz);
            if(Adeltaz < inflray.RadMax || CFG::run::IsTL()) {
                // LP: Changed to use ApplyContribution in order to support
                // incoherent, semi-coherent, and arrivals.
                real cpa = STD::abs(deltaz * (rB - rA))
                    / STD::sqrt(SQ(rB - rA) + SQ(point1.x.y - point0.x.y));
                real ds       = STD::sqrt(SQ(deltaz) - SQ(cpa));
                real sx1      = sint + ds;
                real thet     = STD::atan(cpa / sx1);
                cpx delay     = tau + rayt.y * deltaz;
                real cnst     = inflray.Ratio1 * cn * point1.Amp / STD::sqrt(sx1);
                w             = STD::exp(-a * SQ(thet));
                real phaseInt = point1.Phase + inflray.phase;
                ApplyContribution<CFG, O3D, false>(
                    uAllSources, cnst, w, inflray.omega, delay, phaseInt, RcvrDeclAngle,
                    RL(0.0), 0, inflray.ir, iz, is, inflray, point1, Pos, Beam, eigen,
                    arrinfo);
            }
        }

        inflray.qOld = q;
        ++inflray.ir;
        if(inflray.ir >= Pos->NRr) return false;
    }

    return true;
}

////////////////////////////////////////////////////////////////////////////////
// Main step function
////////////////////////////////////////////////////////////////////////////////

/**
 * LP: Returns whether to continue the ray.
 */
template<typename CFG, bool O3D, bool R3D> HOST_DEVICE inline bool Step_Influence(
    const rayPt<R3D> &point0, const rayPt<R3D> &point1, InfluenceRayInfo<R3D> &inflray,
    int32_t is, cpxf *uAllSources, const BdryType *Bdry, const Origin<O3D, R3D> &org,
    const SSPStructure *ssp, SSPSegState &iSeg, const Position *Pos,
    const BeamStructure<O3D> *Beam, EigenInfo *eigen, const ArrInfo *arrinfo)
{
    if constexpr(R3D) {
        IGNORE_UNUSED(Bdry);
        IGNORE_UNUSED(ssp);
    }

    // See PreRun_Influence, make sure these remain in sync.
    if constexpr(CFG::infl::IsCerveny()) {
        if constexpr(R3D) {
            GlobalLog("Internal error with infl template system");
            bail();
        } else {
            if constexpr(CFG::infl::IsRayCen()) {
                return Step_InfluenceCervenyRayCen<CFG, O3D>(
                    point0, point1, inflray, is, uAllSources, Bdry, Pos, Beam);
            } else {
                return Step_InfluenceCervenyCart<CFG, O3D>(
                    point0, point1, inflray, is, uAllSources, Bdry, org, ssp, iSeg, Pos,
                    Beam);
            }
        }
    } else if constexpr(CFG::infl::IsSGB()) {
        if constexpr(R3D) {
            GlobalLog("Internal error with infl template system");
            bail();
        } else {
            return Step_InfluenceSGB<CFG, O3D>(
                point0, point1, inflray, is, uAllSources, Pos, Beam, eigen, arrinfo);
        }
    } else if constexpr(CFG::infl::IsGeometric()) {
        if constexpr(CFG::infl::IsCartesian()) {
            return Step_InfluenceGeoCart<CFG, O3D, R3D>(
                point0, point1, inflray, is, uAllSources, Pos, Beam, eigen, arrinfo);
        } else {
            return Step_InfluenceGeoRayCen<CFG, O3D, R3D>(
                point0, point1, inflray, is, uAllSources, Pos, Beam, eigen, arrinfo);
        }
    } else {
        static_assert(!sizeof(CFG), "Invalid template in Step_Influence");
    }
}

////////////////////////////////////////////////////////////////////////////////
// Post-processing
////////////////////////////////////////////////////////////////////////////////

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
template<bool O3D, bool R3D> HOST_DEVICE inline void ScalePressure(
    real Dalpha, real Dbeta, real c, cpx epsilon1, cpx epsilon2, float *r, cpxf *u,
    int32_t Ntheta, int32_t NRz, int32_t Nr, real freq, const BeamStructure<O3D> *Beam)
{
    real cnst;
    cpx cnst3d;

    // Compute scale factor for field
    if constexpr(R3D) {
        if(IsCervenyInfl(Beam) && IsCartesianInfl(Beam)) {
            // Cerveny Gaussian beams in Cartesian coordinates
            // epsilon is normally imaginary here, so cnst is complex
            real sqrtc  = STD::sqrt(c);
            real sqrtc3 = CUBE(sqrtc);
            // put this factor into the beam instead?
            cnst3d = STD::sqrt(epsilon1 * epsilon2) * freq * Dbeta * Dalpha / sqrtc3;
            // LP: This is applied before the intensity to pressure conversion.
            for(int32_t itheta = 0; itheta < Ntheta; ++itheta) {
                for(int32_t irz = 0; irz < NRz; ++irz) {
                    for(int32_t ir = 0; ir < Nr; ++ir) {
                        size_t addr = ((size_t)itheta * NRz + irz) * Nr + ir;
                        u[addr] *= Cpx2Cpxf(cnst3d);
                    }
                }
            }
        } else {
            cnst3d = cpx(FL(1.0), FL(0.0)); // LP: set but not used
        }
    } else {
        IGNORE_UNUSED(Dbeta);
        IGNORE_UNUSED(epsilon1);
        IGNORE_UNUSED(epsilon2);
        if(IsCervenyInfl(Beam)) {
            cnst = -Dalpha * STD::sqrt(freq) / c;
        } else {
            cnst = FL(-1.0);
        }
    }

    // For incoherent run, convert intensity to pressure
    if(!IsCoherentRun(Beam)) {
        for(int32_t itheta = 0; itheta < Ntheta; ++itheta) {
            for(int32_t irz = 0; irz < NRz; ++irz) {
                for(int32_t ir = 0; ir < Nr; ++ir) {
                    size_t addr = ((size_t)itheta * NRz + irz) * Nr + ir;
                    u[addr]     = cpxf(STD::sqrt(u[addr].real()), 0.0f);
                }
            }
        }
    }

    if constexpr(!R3D) {
        // scale and/or incorporate cylindrical spreading
        for(int32_t ir = 0; ir < Nr; ++ir) {
            real factor;
            if(IsLineSource(Beam)) {
                const float local_pi = 3.14159265f;
                factor               = FL(-4.0) * STD::sqrt(local_pi) * cnst;
            } else {
                if(r[ir] == 0.0f) {
                    factor = RL(0.0); // avoid /0 at origin, return pressure = 0
                } else {
                    factor = cnst / (real)STD::sqrt(STD::abs(r[ir]));
                }
            }
            for(int32_t irz = 0; irz < NRz; ++irz) { u[irz * Nr + ir] *= (float)factor; }
        }
    }
}

} // namespace bhc
