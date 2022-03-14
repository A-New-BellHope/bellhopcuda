#pragma once
#include "common.hpp"
#include "curves.hpp"

namespace bhc {

constexpr real betaPowerLaw = FL(1.0);

#define SSP_FN_ARGS const vec2 &x, const vec2 &t, cpx &ccpx, vec2 &gradc, \
    real &crr, real &crz, real &czz, real &rho, real freq, \
    const SSPStructure *ssp, int32_t &iSegz, int32_t &iSegr
#define SSP_CALL_ARGS x, t, ccpx, gradc, crr, crz, czz, rho, freq, ssp, iSegz, iSegr
#define SSP_INIT_ARGS vec2 x, const real &fT, \
    LDIFile &ENVFile, std::ostream &PRTFile, std::string FileRoot, \
    SSPStructure *ssp, const AttenInfo *atten, const FreqInfo *freqinfo, HSInfo &RecycledHS
#define SSP_CALL_INIT_ARGS x, fT, ENVFile, PRTFile, FileRoot, ssp, atten, freqinfo, RecycledHS

HOST_DEVICE inline void UpdateDepthSegmentT(const vec2 &x, const vec2 &t,
    const SSPStructure *ssp, int32_t &iSegz)
{
    //LP: Handles edge cases based on which direction the ray is going. If the
    //ray takes a small step in the direction of t, it will remain in the same
    //segment as it is now.
    if(t.y >= RL(0.0)){
        //ssp->z[iSegz] <= x.y < ssp->z[iSegz+1]
        while(x.y < ssp->z[iSegz] && iSegz > 0) --iSegz;
        while(x.y >= ssp->z[iSegz+1] && iSegz < ssp->NPts-2) ++iSegz;
    }else{
        //ssp->z[iSegz] < x.y <= ssp->z[iSegz+1]
        while(x.y > ssp->z[iSegz+1] && iSegz < ssp->NPts-2) ++iSegz;
        while(x.y <= ssp->z[iSegz] && iSegz > 0) --iSegz;
    }
}

HOST_DEVICE inline void UpdateRangeSegmentT(const vec2 &x, const vec2 &t,
    const SSPStructure *ssp, int32_t &iSegr)
{
    //LP: Handles edge cases based on which direction the ray is going. If the
    //ray takes a small step in the direction of t, it will remain in the same
    //segment as it is now.
    if(t.x >= RL(0.0)){
        //ssp->Seg.r[iSegr] <= x.x < ssp->Seg.r[iSegr+1]
        while(x.x < ssp->Seg.r[iSegr] && iSegr > 0) --iSegr;
        while(x.x >= ssp->Seg.r[iSegr+1] && iSegr < ssp->Nr-2) ++iSegr;
    }else{
        //ssp->Seg.r[iSegr] < x.x <= ssp->Seg.r[iSegr+1]
        while(x.x > ssp->Seg.r[iSegr+1] && iSegr < ssp->Nr-2) ++iSegr;
        while(x.x <= ssp->Seg.r[iSegr] && iSegr > 0) --iSegr;
    }
}

HOST_DEVICE inline real LinInterpDensity(const vec2 &x,
    const SSPStructure *ssp, const int32_t &iSegz, real &rho)
{
    real w = (x.y - ssp->z[iSegz]) / (ssp->z[iSegz+1] - ssp->z[iSegz]);
    rho = (RL(1.0) - w) * ssp->rho[iSegz] + w * ssp->rho[iSegz+1];
    return w;
}

/**
 * N2-linear interpolation of SSP data
 */
HOST_DEVICE inline void n2Linear(SSP_FN_ARGS)
{
    IGNORE_UNUSED(freq);
    IGNORE_UNUSED(iSegr);
    
    UpdateDepthSegmentT(x, t, ssp, iSegz);
    real w = LinInterpDensity(x, ssp, iSegz, rho);
    
    ccpx = RL(1.0) / STD::sqrt((RL(1.0) - w) * ssp->n2[iSegz] + w * ssp->n2[iSegz+1]);
    real c = ccpx.real();
    
    gradc = vec2(RL(0.0), RL(-0.5) * c * c * c * ssp->n2z[iSegz].real());
    crr = crz = RL(0.0);
    czz = RL(3.0) * gradc.y * gradc.y / c;
}

/**
 * c-linear interpolation of SSP data
 */
HOST_DEVICE inline void cLinear(SSP_FN_ARGS)
{
    IGNORE_UNUSED(freq);
    IGNORE_UNUSED(iSegr);
    
    UpdateDepthSegmentT(x, t, ssp, iSegz);
    LinInterpDensity(x, ssp, iSegz, rho);
    
    ccpx = ssp->c[iSegz] + (x.y - ssp->z[iSegz]) * ssp->cz[iSegz];
    gradc = vec2(RL(0.0), ssp->cz[iSegz].real());
    crr = crz = czz = RL(0.0);
}

/**
 * This implements the new monotone piecewise cubic Hermite interpolating
 * polynomial (PCHIP) algorithm for the interpolation of the sound speed c.
 */
HOST_DEVICE inline void cPCHIP(SSP_FN_ARGS)
{
    IGNORE_UNUSED(freq);
    IGNORE_UNUSED(iSegr);
    
    UpdateDepthSegmentT(x, t, ssp, iSegz);
    LinInterpDensity(x, ssp, iSegz, rho);
    
    real xt = x.y - ssp->z[iSegz];
    if(STD::abs(xt) > RL(1.0e10)){
        printf("Invalid xt %g\n", xt);
    }
    for(int32_t i=0; i<4; ++i)
        if(STD::abs(ssp->cCoef[i][iSegz]) > RL(1.0e10))
            printf("Invalid ssp->cCoef[%d][%d] = (%g,%g)\n", i, iSegz,
                ssp->cCoef[i][iSegz].real(), ssp->cCoef[i][iSegz].imag());
    
    ccpx = ssp->cCoef[0][iSegz]
        + (ssp->cCoef[1][iSegz]
        + (ssp->cCoef[2][iSegz]
        +  ssp->cCoef[3][iSegz] * xt) * xt) * xt;
    
    gradc = vec2(RL(0.0), (ssp->cCoef[1][iSegz]
              + (RL(2.0) * ssp->cCoef[2][iSegz]
              +  RL(3.0) * ssp->cCoef[3][iSegz] * xt) * xt).real());
    
    crr = crz = RL(0.0);
    czz = (RL(2.0) * ssp->cCoef[2][iSegz] + RL(6.0) * ssp->cCoef[3][iSegz] * xt).real();
}

/**
 * Cubic spline interpolation
 */
HOST_DEVICE inline void cCubic(SSP_FN_ARGS)
{
    IGNORE_UNUSED(freq);
    IGNORE_UNUSED(iSegr);
    
    UpdateDepthSegmentT(x, t, ssp, iSegz);
    LinInterpDensity(x, ssp, iSegz, rho);
    
    real hSpline = x.y - ssp->z[iSegz];
    cpx czcpx, czzcpx;
    
    SplineALL(ssp->cSpline[0][iSegz], ssp->cSpline[1][iSegz], ssp->cSpline[2][iSegz],
        ssp->cSpline[3][iSegz], hSpline, ccpx, czcpx, czzcpx);
    
    // LP: Only for these conversions, BELLHOP uses DBLE() instead of REAL().
    // The manual for DBLE simply says that it converts the argument to double
    // precision and complex is a valid input, but it doesn't say how that
    // conversion is done. Assuming it does real part rather than magnitude.
    gradc = vec2(RL(0.0), czcpx.real());
    crr = crz = RL(0.0);
    czz = czzcpx.real();
}

/**
 * Bilinear quadrilatteral interpolation of SSP data in 2D
 */
HOST_DEVICE inline void Quad(SSP_FN_ARGS)
{
    IGNORE_UNUSED(freq);
    
    real c1, c2, cz1, cz2, cr, cz, s1, s2, delta_r, delta_z;
    
    if(x.x < ssp->Seg.r[0] || x.x > ssp->Seg.r[ssp->Nr-1]){
        printf("Quad: ray is outside the box where the soundspeed is defined\n");
        bail();
    }
    
    UpdateDepthSegmentT(x, t, ssp, iSegz);
    UpdateRangeSegmentT(x, t, ssp, iSegr);
    LinInterpDensity(x, ssp, iSegz, rho);
    if(iSegz >= ssp->Nz - 1 || iSegr >= ssp->Nr - 1){
        printf("iSeg error in Quad: z %d/%d r %d/%d\n",
            iSegz, ssp->Nz, iSegr, ssp->Nr);
        bail();
    }
    
    // for this depth, x.y, get the sound speed at both ends of the segment
    int32_t Nr = ssp->Nr;
    cz1 = ssp->czMat[iSegz*Nr + iSegr  ];
    cz2 = ssp->czMat[iSegz*Nr + iSegr+1];
    
    s2      = x.y - ssp->z[iSegz];
    delta_z = ssp->z[iSegz + 1] - ssp->z[iSegz];
    
    c1 = ssp->cMat[iSegz*Nr + iSegr  ] + s2 * cz1;
    c2 = ssp->cMat[iSegz*Nr + iSegr+1] + s2 * cz2;
    
    // s1 = proportional distance of x.x in range
    delta_r = ssp->Seg.r[iSegr+1] - ssp->Seg.r[iSegr];
    s1 = (x.x - ssp->Seg.r[iSegr]) / delta_r;
    s1 = bhc::min(s1, RL(1.0)); // force piecewise constant extrapolation for points outside the box
    s1 = bhc::max(s1, RL(0.0)); // "
    
    real c = (RL(1.0) - s1) * c1 + s1 * c2;
    
    // interpolate the attenuation !!!! This will use the wrong segment if the ssp in the envil is sampled at different depths
    s2 = s2 / delta_z; // convert to a proportional depth in the layer
    real cimag = ((RL(1.0) - s2) * ssp->c[iSegz] + s2 * ssp->c[iSegz+1]).imag(); // volume attenuation is taken from the single c(z) profile
    
    ccpx = cpx(c, cimag);
    
    cz = (RL(1.0) - s1) * cz1 + s1 * cz2;
    
    cr  = (c2  - c1 ) / delta_r;
    crz = (cz2 - cz1) / delta_r;
    
    gradc = vec2(cr, cz);
    crr = RL(0.0);
    czz = RL(0.0);
}

HOST_DEVICE inline void Analytic(SSP_FN_ARGS)
{
    IGNORE_UNUSED(t);
    IGNORE_UNUSED(freq);
    IGNORE_UNUSED(ssp);
    IGNORE_UNUSED(iSegr);
    
    real c0, cr, cz, DxtDz, xt;
    
    iSegz = 0;
    c0 = FL(1500.0);
    rho = FL(1.0);
    
    const float unk1 = FL(1300.0);
    const float unk2 = FL(0.00737);
    
    // homogeneous halfspace was removed since BELLHOP needs to get gradc just a little below the boundaries, on ray reflection
    
    //if(x.y < 5000.0){
    xt    = FL(2.0) * (x.y - unk1) / unk1;
    real emxt = STD::exp(-xt);
    DxtDz = FL(2.0) / unk1;
    ccpx  = cpx(c0 * (FL(1.0) + unk2 * (xt - FL(1.0) + emxt)), FL(0.0));
    cz    = c0 * unk2 * (FL(1.0) - emxt) * DxtDz;
    czz   = c0 * unk2 * emxt * SQ(DxtDz);
    //}else{
    // Homogeneous half-space
    //xt    = FL(2.0) * (FL(5000.0) - unk1) / unk1;
    //ccpx  = cpx(c0 * (FL(1.0) + unk2 * (xt - FL(1.0) + emxt)), FL(0.0));
    //cz    = FL(0.0);
    //czz   = FL(0.0);
    //}
    
    cr = FL(0.0);
    gradc = vec2(cr, cz);
    crz = FL(0.0);
    crr = FL(0.0);
}

HOST_DEVICE inline void EvaluateSSP(SSP_FN_ARGS)
{
    switch(ssp->Type){
    case 'N': // N2-linear profile option
        n2Linear(SSP_CALL_ARGS); break;
    case 'C': // C-linear profile option
        cLinear (SSP_CALL_ARGS); break;
    case 'P': // monotone PCHIP ACS profile option
        cPCHIP  (SSP_CALL_ARGS); break;
    case 'S': // Cubic spline profile option
        cCubic  (SSP_CALL_ARGS); break;
    case 'Q':
        Quad    (SSP_CALL_ARGS); break;
    /* case 'H':
        // this is called by BELLHOP3D only once, during READIN
        // possibly the logic should be changed to call EvaluateSSP2D or 3D
        x3 = vec3(RL(0.0), RL(0.0), x.y);
        Hexahedral(x3, c, cimag, gradc_3d, cxx, cyy, czz, cxy, cxz, cyz, rho, freq); break; */
    case 'A': // Analytic profile option
        Analytic(SSP_CALL_ARGS); break;
    default:
        printf("EvaluateSSP: Invalid profile option %c\n", ssp->Type);
        bail();
    }
}

HOST_DEVICE inline void EvaluateSSPCOnly(const vec2 &x, const vec2 &t, cpx &ccpx,
    real freq, const SSPStructure *ssp, int32_t &iSegz, int32_t &iSegr)
{
    vec2 gradc;
    real crr, crz, czz, rho;
    EvaluateSSP(SSP_CALL_ARGS);
}

void InitializeSSP(SSP_INIT_ARGS);
 
}
