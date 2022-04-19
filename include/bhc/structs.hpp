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

#ifndef _BHC_INCLUDED_
#error "This file must be included via #include <bhc/bhc.hpp>!"
#endif

#include <type_traits>

namespace bhc {

template<bool THREED> struct TmplVec23 {};
template<> struct TmplVec23<false> { typedef vec2 type; };
template<> struct TmplVec23<true>  { typedef vec3 type; };
template<typename VEC> HOST_DEVICE inline real &DEP(VEC &v);
template<> HOST_DEVICE inline real &DEP(vec2 &v) { return v.y; }
template<> HOST_DEVICE inline real &DEP(vec3 &v) { return v.z; }
template<typename VEC> HOST_DEVICE inline const real &DEP(const VEC &v);
template<> HOST_DEVICE inline const real &DEP(const vec2 &v) { return v.y; }
template<> HOST_DEVICE inline const real &DEP(const vec3 &v) { return v.z; }
template<bool THREED> struct TmplInt12 {};
template<> struct TmplInt12<false> { typedef int32_t type; };
template<> struct TmplInt12<true>  { typedef int2 type; };

////////////////////////////////////////////////////////////////////////////////
//SSP / Boundary
////////////////////////////////////////////////////////////////////////////////

constexpr int32_t MaxN = 100000;
constexpr int32_t MaxSSP = MaxN + 1;

struct rxyz_vector {
    real *r, *x, *y, *z;
};

struct SSPStructure {
    int32_t NPts, Nr, Nx, Ny, Nz;
    real z[MaxSSP], rho[MaxSSP];
    cpx c[MaxSSP], cz[MaxSSP], n2[MaxSSP], n2z[MaxSSP], cSpline[4][MaxSSP];
    cpx cCoef[4][MaxSSP], CSWork[4][MaxSSP]; // for PCHIP coefs.
    real *cMat, *czMat; // LP: No need for separate cMat3 / czMat3 as we don't have to specify the dimension here.
    rxyz_vector Seg;
    char Type;
    char AttenUnit[2];
};

struct SSPOutputs2DExtras {
    real crr, crz;
};
struct SSPOutputs3DExtras {
    real cxx, cyy, cxy, cxz, cyz;
};
template<bool THREED> struct SSPOutputs 
    : std::conditional_t<THREED, SSPOutputs3DExtras, SSPOutputs2DExtras>
{
    cpx ccpx;
    typename TmplVec23<THREED>::type gradc;
    real rho, czz;
};

struct SSPSegState {
    int32_t r, x, y, z;
};

struct HSInfo {
    real alphaR, betaR, alphaI, betaI; // compressional and shear wave speeds/attenuations in user units
    cpx cP, cS; // P-wave, S-wave speeds
    real rho, Depth; // density, depth
    char bc; // Boundary condition type
    char Opt[6];
};

struct BdryPtSmall {
    HSInfo hs;
};

/**
 * LP: There are two boundary structures. This one actually existed in the
 * FORTRAN, and is called Bdry. The other is BdryInfo bdinfo (below).
 */
struct BdryType {
    BdryPtSmall Top, Bot;
};

struct BdryPtFull2DExtras {
    vec2 Nodet; // tangent at the node, if the curvilinear option is used
    real kappa; // curvature of a segement
    real Dx, Dxx, Dss; // first, second derivatives wrt depth; s is along tangent
    HSInfo hs;
};
struct BdryPtFull3DExtras {
    vec3 n1, n2; // (outward-pointing) normals for each of the triangles in a pair, n is selected from those
    vec3 Noden_unscaled;
    real kappa_xx, kappa_xy, kappa_yy;
    real z_xx, z_xy, z_yy;
    real phi_xx, phi_xy, phi_yy;
};
template<bool THREED> struct BdryPtFull 
    : public std::conditional_t<THREED, BdryPtFull3DExtras, BdryPtFull2DExtras>
{
    using VEC = typename TmplVec23<THREED>::type;
    VEC x; // 2D: coordinate for a segment / 3D: coordinate of boundary
    VEC t; // 2D: tangent for a segment / 3D: tangent for a facet
    VEC n; // 2D: outward normal for a segment / 3D: normal for a facet (outward pointing)
    VEC Noden; // normal at the node, if the curvilinear option is used
    real Len; // 2D: length of a segment / 3D: length of tangent (temporary variable to normalize tangent)
};

template<bool THREED> struct BdryInfoTopBot {
    using IORI2 = typename TmplInt12<THREED>::type;
    IORI2 NPts;
    char type[2]; // In 3D, only first char is used
    BdryPtFull<THREED> *bd; // 2D: 1D array / 3D: 2D array
};
/**
 * LP: There are two boundary structures. This one represents data in the
 * FORTRAN which was not within any structure, and is called bdinfo.
 */
template<bool THREED> struct BdryInfo {
    BdryInfoTopBot<THREED> top, bot;
};

struct BdryLimits { real min, max; };
struct BdryLimits2 { BdryLimits x, y; };
template<bool THREED> struct TmplBdryLimits12 {};
template<> struct TmplBdryLimits12<false> { typedef BdryLimits type; };
template<> struct TmplBdryLimits12<true>  { typedef BdryLimits2 type; };
template<bool THREED> struct BdryStateTopBot {
    using IORI2 = typename TmplInt12<THREED>::type;
    using VEC = typename TmplVec23<THREED>::type;
    using LIMITS = typename TmplBdryLimits12<THREED>::type;
    IORI2 Iseg; // indices that point to the current active segment
    LIMITS lSeg; // LP: limits of current segment
    VEC x, n; // only explicitly written in 3D, but effectively present in 2D
};
/**
 * LP: Variables holding current state of active segment(s) ray is in.
 */
template<bool THREED> struct BdryState {
    BdryStateTopBot<THREED> top, bot;
};

////////////////////////////////////////////////////////////////////////////////
//Reflections
////////////////////////////////////////////////////////////////////////////////

struct ReflectionCoef {
    real theta, r, phi;
};

struct ReflectionInfo {
    int32_t NBotPts, NTopPts;
    ReflectionCoef *RBot, *RTop;
};

////////////////////////////////////////////////////////////////////////////////
//Attenuation (volume absorption / scattering)
////////////////////////////////////////////////////////////////////////////////

struct bioStructure {
    real z1, z2, f0, q, a0;
};

constexpr int32_t MaxBioLayers = 200;

struct AttenInfo {
    int32_t NBioLayers;
    bioStructure bio[MaxBioLayers];
    real t, Salinity, pH, z_bar, fg; // Francois-Garrison volume attenuation; temperature, salinity, ...
};

////////////////////////////////////////////////////////////////////////////////
//Source / receiver positions
////////////////////////////////////////////////////////////////////////////////

struct Position {
    int32_t NSx, NSy, NSz, NRz, NRr, Ntheta; // number of x, y, z, r, theta coordinates
    int32_t NRz_per_range;
    real Delta_r, Delta_theta;
    int32_t *iSz, *iRz;
    // LP: These are really floats, not reals.
    float *Sx, *Sy, *Sz; // Source x, y, z coordinates
    float *Rr, *Rz, *ws, *wr; // Receiver r, z coordinates and weights for interpolation
    float *theta; // Receiver bearings
};

////////////////////////////////////////////////////////////////////////////////
//Source angles
////////////////////////////////////////////////////////////////////////////////

struct AnglesStructure {
    int32_t Nalpha, Nbeta, iSingle_alpha, iSingle_beta;
    real Dalpha, Dbeta; // angular spacing between beams
    real *alpha; // LP: elevation angles
    real *beta; // LP: azimuth angles
};

////////////////////////////////////////////////////////////////////////////////
//Source frequencies
////////////////////////////////////////////////////////////////////////////////

struct FreqInfo {
    real freq0; // Nominal or carrier frequency
    int32_t Nfreq; // number of frequencies
    real *freqVec; // frequency vector for braoaband runs
};

////////////////////////////////////////////////////////////////////////////////
//Beams
////////////////////////////////////////////////////////////////////////////////

struct rxyz {
    real r, x, y, z;
};

/**
 * LP: Like boundaries, there are two beam structures. This one is (mostly)
 * in the FORTRAN, and is called Beam.
 */
struct BeamStructure {
    //LP: NSteps moved out of this struct as it's a property of a single beam.
    int32_t NBeams, Nimage, iBeamWindow;
    real deltas, epsMultiplier, rLoop;
    char Component;
    char Type[4];
    char RunType[7];
    rxyz Box;
};

/**
 * LP: TODO: rename to something with SBP, this is not the same kind of beam
 * as BeamStructure.
 */
struct BeamInfo {
    int32_t NSBPPts;
    real *SrcBmPat;
    char SBPFlag;
};

////////////////////////////////////////////////////////////////////////////////
//Eigenrays
////////////////////////////////////////////////////////////////////////////////

struct EigenHit {
    // LP: Receiver this hit pertains to
    int32_t ir, iz;
    // LP: Identifying info to re-trace this ray
    int32_t isrc, ialpha;
    // LP: Number of steps until the ray hit the receiver
    int32_t is;
};

struct EigenInfo {
    uint32_t neigen;
    uint32_t memsize;
    EigenHit *hits;
};

////////////////////////////////////////////////////////////////////////////////
//Arrivals
////////////////////////////////////////////////////////////////////////////////

struct Arrival {
    int32_t NTopBnc, NBotBnc;
    float SrcDeclAngle, SrcAzimAngle, RcvrDeclAngle, RcvrAzimAngle, a, Phase;
    cpxf delay;
};

/**
 * LP: Arrival setup and results.
 */
struct ArrInfo {
    Arrival *Arr;
    int32_t *NArr;
    int32_t MaxNArr;
    bool singlethread;
};

////////////////////////////////////////////////////////////////////////////////
//Rays/beams
////////////////////////////////////////////////////////////////////////////////

struct ray2DPt {
    int32_t NumTopBnc, NumBotBnc;
    ///ray coordinate, (r,z)
    vec2 x;
    ///scaled tangent to the ray (previously (rho, zeta))
    vec2 t;
    vec2 p, q;
    ///c * t would be the unit tangent
    real c;
    real Amp, Phase;
    cpx tau;
};

struct RayResult {
    ray2DPt *ray2D;
    real SrcDeclAngle;
    int32_t Nsteps;
};

struct RayInfo {
    ray2DPt *raymem;
    uint32_t NPoints;
    uint32_t MaxPoints;
    RayResult *results;
    int32_t NRays;
};

////////////////////////////////////////////////////////////////////////////////
//Meta-structures
////////////////////////////////////////////////////////////////////////////////

struct bhcParams {
    char Title[80]; // Size determined by WriteHeader for TL
    real fT;
    BdryType *Bdry;
    BdryInfo<false> *bdinfo;
    ReflectionInfo *refl;
    SSPStructure *ssp;
    AttenInfo *atten;
    Position *Pos;
    AnglesStructure *Angles;
    FreqInfo *freqinfo;
    BeamStructure *Beam;
    BeamInfo *beaminfo;
    ///Pointer to internal data structure for program (non-marine-related) state.
    void *internal;
};

struct bhcOutputs {
    RayInfo *rayinfo;
    cpxf *uAllSources;
    EigenInfo *eigen;
    ArrInfo *arrinfo;
};

}
