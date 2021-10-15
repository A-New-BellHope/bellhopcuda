#pragma once
#include "common.hpp"

struct rxyz {
    real r, x, y, z;
};

struct rxyz_vector {
    real *r, *x, *y, *z;
};

struct BeamStructure {
    int32_t NBeams, Nimage, Nsteps, iBeamWindow;
    real deltas, epsMultiplier = RC(1.0), rLoop;
    char Component;
    char Type[4] = "G S ";
    char RunType[7];
    rxyz Box;
};

constexpr int32_t MaxN = 100000;
constexpr int32_t MaxSSP = MaxN + 1;

struct SSPStructure {
    int32_t NPts, Nr, Nx, Ny, Nz;
    real z[MaxSSP], rho[MaxSSP];
    cpx c[MaxSSP], cz[MaxSSP], n2[MaxSSP], n2z[MaxSSP], cSpline[4][MaxSSP];
    cpx cCoef[4][MaxSSP], CSWork[4][MaxSSP]; // for PCHIP coefs.
    real *cMat, *czMat, *cMat3, *czMat3;
    rxyz_vector Seg;
    char Type;
    char AttenUnit[2];
};

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

struct HSInfo {
    cpx alpha, beta; // compressional and shear wave speeds/attenuations in user units
    cpx cP, cS; // P-wave, S-wave speeds
    real rho, Depth; // density, depth
    char bc; // Boundary condition type
    char[6] Opt;
};

struct BdryPtSmall {
    HSInfo hs;
};

struct BdryPtFull {
    vec2 x, t, n; // coordinate, tangent, and outward normal for a segment
    vec2 Nodet, Noden; // tangent and normal at the node, if the curvilinear option is used
    real Len, Kappa; // length and curvature of a segement
    real Dx, Dxx, Dss; // first, second derivatives wrt depth; s is along tangent
    HSInfo hs;
};

struct BdryType {
    BdryPtSmall Top, Bot;
};

struct ReflectionCoef {
    real theta, r, phi;
};

struct AnglesStructure {
    int32_t Nalpha = 0, Nbeta = 1, iSingle_alpha = -1, iSingle_beta = -1;
    real Dalpha, Dbeta;
    real *alpha, *beta;
};

struct bioStructure {
    real Z1, Z2, f0, q, a0;
};

constexpr int32_t MaxBioLayers = 200;

struct AttenInfo {
    int32_t NBioLayers;
    bioStructure bio[MaxBioLayers];
    real t = RC(20.0), Salinity = RC(35.0), pH = RC(8.0), z_bar = RC(0.0), fg;
};

struct Position {
    int32_t NSx = 1, NSy = 1, NSz, NRz, NRr, Ntheta; // number of x, y, z, r, theta coordinates
    real Delta_r, Delta_theta;
    int32_t *iSz, *iRz;
    real *Sx, *Sy, *Sz; // Source x, y, z coordinates
    real *Rr, *Rz, *ws, *wr; // Receiver r, z coordinates and weights for interpolation
    real *theta; // Receiver bearings
};
