/*
bellhopcxx / bellhopcuda - C++/CUDA port of BELLHOP(3D) underwater acoustics simulator
Copyright (C) 2021-2023 The Regents of the University of California
Marine Physical Lab at Scripps Oceanography, c/o Jules Jaffe, jjaffe@ucsd.edu
Based on BELLHOP / BELLHOP3D, which is Copyright (C) 1983-2022 Michael B. Porter

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

/*
O3D: ocean (SSP, boundaries, etc.) is 3D.
R3D: rays (beams) are 3D. In 2D-3D mode, S3D == true && R3D == false.
X3D: generic 3D boolean value for templated types.
*/

template<bool X3D> struct TmplVec23 {};
template<> struct TmplVec23<false> {
    typedef vec2 type;
};
template<> struct TmplVec23<true> {
    typedef vec3 type;
};
template<bool X3D> using VEC23 = typename TmplVec23<X3D>::type;

template<bool X3D> struct TmplInt12 {};
template<> struct TmplInt12<false> {
    typedef int32_t type;
};
template<> struct TmplInt12<true> {
    typedef int2 type;
};
template<bool X3D> using IORI2 = typename TmplInt12<X3D>::type;

template<bool X3D> struct TmplVec2Mat2 {};
template<> struct TmplVec2Mat2<false> {
    typedef vec2 type;
};
template<> struct TmplVec2Mat2<true> {
    typedef mat2x2 type;
};
template<bool X3D> using V2M2 = typename TmplVec2Mat2<X3D>::type;

////////////////////////////////////////////////////////////////////////////////
// SSP / Boundary
////////////////////////////////////////////////////////////////////////////////

constexpr int32_t MaxN   = 100000;
constexpr int32_t MaxSSP = MaxN + 1;

struct rxyz_vector {
    real *r, *x, *y, *z;
};

struct SSPStructure {
    // LP: Start with complex values for alignment reasons.
    cpx c[MaxSSP], cz[MaxSSP], n2[MaxSSP], n2z[MaxSSP], cSpline[4][MaxSSP];
    cpx cCoef[4][MaxSSP], CSWork[4][MaxSSP]; // for PCHIP coefs.
    real *cMat, *czMat; // LP: No need for separate cMat3 / czMat3 as we don't have to
                        // specify the dimension here.
    rxyz_vector Seg;
    real z[MaxSSP], rho[MaxSSP];
    real alphaR[MaxSSP], alphaI[MaxSSP];
    // LP: Not actually used, but echoed, so with new system need to store them
    real betaR[MaxSSP], betaI[MaxSSP];

    int32_t NPts, Nr, Nx, Ny, Nz;
    char Type;
    char AttenUnit[2];
    bool rangeInKm; // Ranges (R, X, Y) specified in km, will be automatically converted
                    // to meters
    bool dirty;     // reset and update derived params
};

// SSPOutputs is templated as O3D during computation, but R3D when returned,
// except in RayStartNominalSSP and everywhere which uses those results.
template<bool X3D> struct SSPOutputsExtras {};
template<> struct SSPOutputsExtras<false> {
    real crr, crz;
};
template<> struct SSPOutputsExtras<true> {
    real cxx, cyy, cxy, cxz, cyz;
};
template<bool X3D> struct SSPOutputs : public SSPOutputsExtras<X3D> {
    VEC23<X3D> gradc;
    cpx ccpx;
    real rho, czz;
};

struct SSPSegState {
    int32_t r, x, y, z;
};

struct HSInfo {
    real alphaR, betaR, alphaI, betaI; // compressional and shear wave speeds/attenuations
                                       // in user units
    cpx cP, cS;                        // P-wave, S-wave speeds
    real rho, Depth;                   // density, depth
    char bc;                           // Boundary condition type
    char Opt[6];
};

struct HSExtra {
    real zTemp, Mz;
    real Sigma; // Read and echoed but never used
};

struct BdryPtSmall {
    HSInfo hs;
    HSExtra hsx;
};

/**
 * LP: There are three boundary structures. This one actually existed in the
 * FORTRAN, and is called Bdry. This holds the boundary properties for the
 * current segment.
 */
struct BdryType {
    BdryPtSmall Top, Bot;
};

template<bool O3D> struct ReflCurvature {};
template<> struct ReflCurvature<false> {
    real kappa; // curvature of a segment
};
template<> struct ReflCurvature<true> {
    real z_xx, z_xy, z_yy, kappa_xx, kappa_xy, kappa_yy;
};

template<bool O3D> struct BdryPtFullExtras {};
template<> struct BdryPtFullExtras<false> {
    vec2 Nodet;        // tangent at the node, if the curvilinear option is used
    real Dx, Dxx, Dss; // first, second derivatives wrt depth; s is along tangent
    HSInfo hs;
};
template<> struct BdryPtFullExtras<true> {
    vec3 n1, n2; // (outward-pointing) normals for each of the triangles in a pair, n is
                 // selected from those
    vec3 Noden_unscaled;
    real phi_xx, phi_xy, phi_yy;
};
template<bool O3D> struct BdryPtFull : public BdryPtFullExtras<O3D>,
                                       public ReflCurvature<O3D> {
    VEC23<O3D> x; // 2D: coordinate for a segment / 3D: coordinate of boundary
    VEC23<O3D> t; // 2D: tangent for a segment / 3D: tangent for a facet
    VEC23<O3D> n; // 2D: outward normal for a segment / 3D: normal for a facet (outward
                  // pointing)
    VEC23<O3D> Noden; // normal at the node, if the curvilinear option is used
    real Len; // 2D: length of a segment / 3D: length of tangent (temporary variable to
              // normalize tangent)
};
template<bool O3D> constexpr int32_t BdryStride = sizeof(BdryPtFull<O3D>) / sizeof(real);

template<bool O3D> struct BdryInfoTopBot {
    IORI2<O3D> NPts;
    char type[2];        // In 3D, only first char is used
    bool dirty;          // Set to indicate that derived values need updating
    bool rangeInKm;      // R, X, Y values in km; automatically converted to meters
    BdryPtFull<O3D> *bd; // 2D: 1D array / 3D: 2D array
};
/**
 * LP: There are three boundary structures. This one represents static/global
 * data in the FORTRAN which was not within any structure, and is called bdinfo.
 */
template<bool O3D> struct BdryInfo {
    BdryInfoTopBot<O3D> top, bot;
};

struct BdryTriDiagState {
    bool side;
    bool onEdge;
    bool justSteppedTo;
    bool outgoingSide;
};

struct BdryLimits {
    real min, max;
};
struct BdryLimits2 {
    BdryLimits x, y;
};
template<bool O3D> struct TmplBdryLimits12 {};
template<> struct TmplBdryLimits12<false> {
    typedef BdryLimits type;
};
template<> struct TmplBdryLimits12<true> {
    typedef BdryLimits2 type;
};
template<bool O3D> struct BdryStateTopBotExtras {};
template<> struct BdryStateTopBotExtras<true> {
    BdryTriDiagState td;
};
template<bool O3D> struct BdryStateTopBot : public BdryStateTopBotExtras<O3D> {
    IORI2<O3D> Iseg; // indices that point to the current active segment
    typename TmplBdryLimits12<O3D>::type lSeg; // LP: limits of current segment
    VEC23<O3D> x, n; // only explicitly written in 3D, but effectively present in 2D
    VEC23<O3D> xmid; // coordinates of center of active rectangle (3D) / segment (2D)
    // because corners may be at big number and mess up floating point precision
};
/**
 * LP: There are three boundary structures. This one holds variables describing
 * where the ray currently is in terms of boundary segments, and is called bds.
 */
template<bool O3D> struct BdryState {
    BdryStateTopBot<O3D> top, bot;
};

/**
 * LP: In Nx2D mode, describes the position of the 2D ray space relative to the
 * 3D ocean. Unused (empty struct) in 2D and full 3D mode.
 */
template<bool O3D, bool R3D> struct Origin {};
template<> struct Origin<true, false> {
    vec3 xs;
    vec2 tradial;
};

////////////////////////////////////////////////////////////////////////////////
// Reflections
////////////////////////////////////////////////////////////////////////////////

struct ReflectionCoef {
    real theta, r, phi;
};

struct ReflectionInfoTopBot {
    int32_t NPts;
    bool inDegrees; // Angles in degrees, converted to radians at preprocess
    ReflectionCoef *r;
};
struct ReflectionInfo {
    ReflectionInfoTopBot bot, top;
};

////////////////////////////////////////////////////////////////////////////////
// Attenuation (volume absorption / scattering)
////////////////////////////////////////////////////////////////////////////////

struct bioStructure {
    real z1, z2, f0, q, a0;
};

constexpr int32_t MaxBioLayers = 200;

struct AttenInfo {
    int32_t NBioLayers;
    bioStructure bio[MaxBioLayers];
    real t, Salinity, pH, z_bar, fg; // Francois-Garrison volume attenuation; temperature,
                                     // salinity, ...
    int32_t NMedia;
};

////////////////////////////////////////////////////////////////////////////////
// Source / receiver positions
////////////////////////////////////////////////////////////////////////////////

struct Position {
    int32_t NSx, NSy, NSz, NRz, NRr, Ntheta; // number of x, y, z, r, theta coordinates
    int32_t NRz_per_range;
    bool SxSyInKm, RrInKm; // Values in km, converted to meters in preprocess
    /// Whether a duplicate angle (e.g. 360.0 when 0.0 also exists) was removed
    /// while reading the environment file. If manually setting theta, set this
    /// to false and do not include any duplicate angles.
    bool thetaDuplRemoved;
    real Delta_r, Delta_theta;
    // int32_t *iSz, *iRz; // LP: Not used.
    // LP: These are really floats, not reals.
    float *Sx, *Sy, *Sz; // Source x, y, z coordinates
    float *Rr, *Rz;      // Receiver r, z coordinates
    // float *ws, *wr; // weights for interpolation LP: Not used.
    float *theta; // Receiver bearings
    vec2 *t_rcvr; // Receiver directions (cos(theta), sin(theta))
};

////////////////////////////////////////////////////////////////////////////////
// Source angles
////////////////////////////////////////////////////////////////////////////////

struct AngleInfo {
    int32_t n, iSingle;
    real d; // angular spacing between beams
    real *angles;
    bool inDegrees; // Angles in degrees, converted to radians at preprocess
};

struct AnglesStructure {
    AngleInfo alpha; // LP: elevation angles
    AngleInfo beta;  // LP: azimuth angles
};

////////////////////////////////////////////////////////////////////////////////
// Source frequencies
////////////////////////////////////////////////////////////////////////////////

struct FreqInfo {
    real freq0;    // Nominal or carrier frequency
    int32_t Nfreq; // number of frequencies
    real *freqVec; // frequency vector for broadband runs
};

////////////////////////////////////////////////////////////////////////////////
// Beams
////////////////////////////////////////////////////////////////////////////////

/**
 * LP: Like boundaries, there are two beam structures. This one is (mostly)
 * in the FORTRAN, and is called Beam.
 */
template<bool O3D> struct BeamStructure {
    // LP: Nsteps moved out of this struct as it's a property of a single beam.
    int32_t NBeams, Nimage, iBeamWindow;
    char Component;
    char Type[4];
    char RunType[7];
    bool rangeInKm;  // Box R, X, Y specified in km, converted to meters in preprocess
    bool autoDeltas; // stores whether deltas was automatically computed, for echo
    real deltas, epsMultiplier, rLoop;
    VEC23<O3D> Box;
};

/**
 * LP: Source Beam Pattern info.
 */
struct SBPInfo {
    char SBPFlag;
    bool SBPIndB; // SrcBmPat values in dB, auto converted to linear in preprocess
    int32_t NSBPPts;
    real *SrcBmPat;
};

////////////////////////////////////////////////////////////////////////////////
// Eigenrays
////////////////////////////////////////////////////////////////////////////////

struct EigenHit {
    // LP: Number of steps until the ray hit the receiver
    int32_t is;
    // LP: Receiver this hit pertains to
    int32_t itheta, ir, iz;
    // LP: Identifying info to re-trace this ray
    int32_t isx, isy, isz, ialpha, ibeta;
};

struct EigenInfo {
    int32_t neigen;
    int32_t memsize;
    EigenHit *hits;
};

////////////////////////////////////////////////////////////////////////////////
// Arrivals
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
    int32_t *MaxNPerSource;
    int32_t MaxNArr;
    bool AllowMerging;
};

////////////////////////////////////////////////////////////////////////////////
// Rays/beams
////////////////////////////////////////////////////////////////////////////////

template<bool R3D> struct rayPtExtras {};
template<> struct rayPtExtras<false> {
    vec2 p, q;
};
template<> struct rayPtExtras<true> {
    // LP: p and q don't actually affect the ray trajectory or results; they are
    // used to update the corresponding variables for the next step, but that's
    // never actually used for anything else. There is commented out code which
    // uses them, as well as additional commented-out variables in this struct,
    // for 3D Cerveny beams.
    mat2x2 p, q; // LP: The ROWS of p are p_tilde and p_hat; same for q.
    // real DetQ; LP: Precomputed at the start of Influence functions, only used
    // for detecting phase inversions. Changed to compute when needed.
    real phi;
    // MSVC warns about rayPt<true> not being aligned for the cpx member
    real dummy;
};
template<bool R3D> struct rayPt : public rayPtExtras<R3D> {
    int32_t NumTopBnc, NumBotBnc;
    /// ray coordinate, (r,z)
    VEC23<R3D> x;
    /// scaled tangent to the ray (previously (rho, zeta))
    VEC23<R3D> t;
    /// c * t would be the unit tangent
    real c;
    real Amp, Phase;
    cpx tau;
};

template<bool R3D> struct StepPartials {};
template<> struct StepPartials<false> {
    real cnn_csq;
};
template<> struct StepPartials<true> {
    real cnn, cmn, cmm;
};

template<bool O3D, bool R3D> struct RayResult {
    rayPt<R3D> *ray;
    Origin<O3D, R3D> org;
    real SrcDeclAngle;
    int32_t Nsteps;
};

template<bool O3D, bool R3D> struct RayInfo {
    RayResult<O3D, R3D> *results;
    rayPt<R3D> *RayMem;
    rayPt<R3D> *WorkRayMem;
    size_t RayMemCapacity;
    size_t RayMemPoints;
    int32_t MaxPointsPerRay;
    int32_t NRays;
    bool isCopyMode;
};

struct RayInitInfo {
    int32_t isx, isy, isz, ialpha, ibeta;
    real alpha, beta;
    real SrcDeclAngle, SrcAzimAngle;
};

////////////////////////////////////////////////////////////////////////////////
// Influence / transmission loss
////////////////////////////////////////////////////////////////////////////////

template<bool R3D> struct InfluenceRayInfo {
    // LP: Constants.
    RayInitInfo init;
    real Dalpha, Dbeta;     // angular spacing
    real c0;                // LP: c at start of ray
    cpx epsilon1, epsilon2; // beam constant
    VEC23<R3D> xs;          // source
    real freq0, omega;
    real RadMax;
    real BeamWindow;
    int32_t iBeamWindow2;
    real Ratio1; // scale factor (point source vs. line source)
    real rcp_q0, rcp_qhat0;
    // LP: Variables carried over between iterations.
    real phase;
    real qOld;               // LP: Det_QOld in 3D
    VEC23<R3D> rayn1, rayn2; // LP: rayn1 was rn, zn in 2D
    VEC23<R3D> x;
    cpx gamma;
    bool lastValid;
    int32_t kmah;
    int32_t ir;
};

////////////////////////////////////////////////////////////////////////////////
// Meta-structures
////////////////////////////////////////////////////////////////////////////////

struct bhcInit {
    /// Number of worker threads to run. -1 means "all logical cores".
    int32_t numThreads = -1;
    /// Maximum amount of memory (in bytes) this instance should use.
    size_t maxMemory = 4ull * 1024ull * 1024ull * 1024ull; // 4 GiB
    /// If there is not enough memory to hold the requested number of rays
    /// where each is maximum length, whether to solve this by reducing the
    /// maximum length (false), or by using copy mode (true). Copy mode can fit
    /// more ray data in memory but is slower. This only affects ray and
    /// eigenray runs (no effect on TL or arrivals).
    bool useRayCopyMode = false;
    /// Index of the GPU to use (ignored if not in CUDA mode). This is the order
    /// the GPUs are enumerated in CUDA, usually with the most powerful GPU
    /// as index 0.
    int gpuIndex = 0;
    /**
     * If not null: Relative path to environment file, without the .env
     * extension. E.g. path/to/MunkB_ray_rot (where path/to/MunkB_ray_rot.env
     * and also path/to/MunkB_ray_rot.ssp, path/to/MunkB_ray_rot.bty, etc.
     * exist). The string passed here will be copied internally and may be
     * deleted by the caller after bhc::setup() returns.

     * If null: Sets up the environment with some very basic defaults.
     * prtCallback must not be null, because if there is no environment file, a
     * real *.prt file cannot be created. bhc::writeout() also cannot be used.
     */
    const char *FileRoot = nullptr;
    /**
     * prtCallback, outputCallback: There are two different types of output
     * messages which can be produced by BELLHOP(3D) and therefore bellhopcxx /
     * bellhopcuda: outputs which form the PRTFile (print file *.prt), and
     * outputs which are printed in the terminal (usually fatal error messages).
     * Two callbacks are provided here to receive these messages; you may handle
     * them or pass nullptr for either/both. If prtCallback is nullptr, a *.prt
     * file will be created with the PRTFile outputs. If outputCallback is
     * nullptr, these messages will be printed to standard output.
     *
     * If you are using multiple instances (multiple calls to setup with
     * different params), any callback you pass for either of these must be
     * thread-safe, as it may be called by multiple threads in parallel.
     * Furthermore, if you are using multiple instances and you set prtCallback
     * to nullptr so it writes PRTFiles, each instance must use a different
     * FileRoot or there will be issues with the multiple instances trying to
     * write to the same PRTFile.
     *
     * The memory pointed to by message will be freed immediately after the
     * callback function returns, so your callback function must copy the string
     * rather than storing the pointer.
     */
    void (*prtCallback)(const char *message) = nullptr;
    /// See documentation for prtCallback above.
    void (*outputCallback)(const char *message) = nullptr;
};

template<bool O3D> struct bhcParams {
    char Title[80]; // Size determined by WriteHeader for TL
    real fT;
    BdryType *Bdry;
    BdryInfo<O3D> *bdinfo;
    ReflectionInfo *refl;
    SSPStructure *ssp;
    AttenInfo *atten;
    Position *Pos;
    AnglesStructure *Angles;
    FreqInfo *freqinfo;
    BeamStructure<O3D> *Beam;
    SBPInfo *sbp;
    /// Pointer to internal data structure for program (non-marine-related) state.
    void *internal;
};

template<bool O3D, bool R3D> struct bhcOutputs {
    RayInfo<O3D, R3D> *rayinfo;
    cpxf *uAllSources;
    EigenInfo *eigen;
    ArrInfo *arrinfo;
};

} // namespace bhc
