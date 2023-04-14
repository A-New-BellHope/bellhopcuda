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
#include "common_run.hpp"
#include "curves.hpp"

namespace bhc {

#define SSP_2D_FN_ARGS \
    const vec2 &x, const vec2 &t, SSPOutputs<false> &o, const SSPStructure *ssp, \
        SSPSegState &iSeg, [[maybe_unused]] ErrState *errState
#define SSP_3D_FN_ARGS \
    const vec3 &x, const vec3 &t, SSPOutputs<true> &o, const SSPStructure *ssp, \
        SSPSegState &iSeg, ErrState *errState
#define SSP_TEMPL_FN_ARGS \
    const VEC23<O3D> &x, [[maybe_unused]] const VEC23<O3D> &t, SSPOutputs<O3D> &o, \
        [[maybe_unused]] const SSPStructure *ssp, SSPSegState &iSeg, \
        [[maybe_unused]] ErrState *errState

HOST_DEVICE inline void UpdateSSPSegment(
    real x, real t, const real *array, int32_t n, int32_t &iSeg)
{
    // LP: Handles edge cases based on which direction the ray is going. If the
    // ray takes a small step in the direction of t, it will remain in the same
    // segment as it is now.
    if(t >= RL(0.0)) {
        // array[iSeg] <= x < array[iSeg+1]
        while(iSeg > 0 && x < array[iSeg]) --iSeg;
        while(iSeg < n - 2 && x >= array[iSeg + 1]) ++iSeg;
    } else {
        // array[iSeg] < x <= array[iSeg+1]
        while(iSeg < n - 2 && x > array[iSeg + 1]) ++iSeg;
        while(iSeg > 0 && x <= array[iSeg]) --iSeg;
    }
}

HOST_DEVICE inline real LinInterpDensity(
    real z, const SSPStructure *ssp, const SSPSegState &iSeg, real &rho)
{
    real w = (z - ssp->z[iSeg.z]) / (ssp->z[iSeg.z + 1] - ssp->z[iSeg.z]);
    rho    = (RL(1.0) - w) * ssp->rho[iSeg.z] + w * ssp->rho[iSeg.z + 1];
    return w;
}

/**
 * N2-linear interpolation of SSP data
 */
HOST_DEVICE inline void n2Linear(SSP_2D_FN_ARGS)
{
    UpdateSSPSegment(x.y, t.y, ssp->z, ssp->NPts, iSeg.z);
    real w = LinInterpDensity(x.y, ssp, iSeg, o.rho);

    o.ccpx = RL(1.0)
        / STD::sqrt((RL(1.0) - w) * ssp->n2[iSeg.z] + w * ssp->n2[iSeg.z + 1]);
    real c = o.ccpx.real();

    o.gradc = vec2(RL(0.0), RL(-0.5) * CUBE(c) * ssp->n2z[iSeg.z].real());
    o.crr = o.crz = RL(0.0);
    o.czz         = RL(3.0) * o.gradc.y * o.gradc.y / c;
}

/**
 * c-linear interpolation of SSP data
 */
HOST_DEVICE inline void cLinear(SSP_2D_FN_ARGS)
{
    UpdateSSPSegment(x.y, t.y, ssp->z, ssp->NPts, iSeg.z);
    LinInterpDensity(x.y, ssp, iSeg, o.rho);

    o.ccpx  = ssp->c[iSeg.z] + (x.y - ssp->z[iSeg.z]) * ssp->cz[iSeg.z];
    o.gradc = vec2(RL(0.0), ssp->cz[iSeg.z].real());
    o.crr = o.crz = o.czz = RL(0.0);
}

/**
 * This implements the new monotone piecewise cubic Hermite interpolating
 * polynomial (PCHIP) algorithm for the interpolation of the sound speed c.
 */
HOST_DEVICE inline void cPCHIP(SSP_2D_FN_ARGS)
{
    UpdateSSPSegment(x.y, t.y, ssp->z, ssp->NPts, iSeg.z);
    LinInterpDensity(x.y, ssp, iSeg, o.rho);

    real xt = x.y - ssp->z[iSeg.z];
    if(STD::abs(xt) > RL(1.0e10)) {
        RunWarning(errState, BHC_WARN_CPCHIP_INVALIDXT);
        // printf("Invalid xt %g\n", xt);
    }
    for(int32_t i = 0; i < 4; ++i) {
        if(STD::abs(ssp->cCoef[i][iSeg.z]) > RL(1.0e10)) {
            RunWarning(errState, BHC_WARN_CPCHIP_INVALIDCCOEF);
            // printf(
            //     "Invalid ssp->cCoef[%d][%d] = (%g,%g)\n", i, iSeg.z,
            //     ssp->cCoef[i][iSeg.z].real(), ssp->cCoef[i][iSeg.z].imag());
        }
    }

    o.ccpx = ssp->cCoef[0][iSeg.z]
        + (ssp->cCoef[1][iSeg.z]
           + (ssp->cCoef[2][iSeg.z] + ssp->cCoef[3][iSeg.z] * xt) * xt)
            * xt;

    o.gradc = vec2(
        RL(0.0),
        (ssp->cCoef[1][iSeg.z]
         + (RL(2.0) * ssp->cCoef[2][iSeg.z] + RL(3.0) * ssp->cCoef[3][iSeg.z] * xt) * xt)
            .real());

    o.crr = o.crz = RL(0.0);
    o.czz
        = (RL(2.0) * ssp->cCoef[2][iSeg.z] + RL(6.0) * ssp->cCoef[3][iSeg.z] * xt).real();
}

/**
 * Cubic spline interpolation
 */
HOST_DEVICE inline void cCubic(SSP_2D_FN_ARGS)
{
    UpdateSSPSegment(x.y, t.y, ssp->z, ssp->NPts, iSeg.z);
    LinInterpDensity(x.y, ssp, iSeg, o.rho);

    real hSpline = x.y - ssp->z[iSeg.z];
    cpx czcpx, czzcpx;

    SplineALL(
        ssp->cSpline[0][iSeg.z], ssp->cSpline[1][iSeg.z], ssp->cSpline[2][iSeg.z],
        ssp->cSpline[3][iSeg.z], hSpline, o.ccpx, czcpx, czzcpx);

    // LP: Only for these conversions, BELLHOP uses DBLE() instead of REAL().
    // The manual for DBLE simply says that it converts the argument to double
    // precision and complex is a valid input, but it doesn't say how that
    // conversion is done. Assuming it does real part rather than magnitude.
    o.gradc = vec2(RL(0.0), czcpx.real());
    o.crr = o.crz = RL(0.0);
    o.czz         = czzcpx.real();
}

/**
 * Bilinear quadrilateral interpolation of SSP data in 2D
 */
HOST_DEVICE inline void Quad(SSP_2D_FN_ARGS)
{
    real c1, c2, cz1, cz2, cr, cz, s1, s2, delta_r, delta_z;

    if(x.x < ssp->Seg.r[0] || x.x > ssp->Seg.r[ssp->Nr - 1]) {
        RunError(errState, BHC_ERR_OUTSIDE_SSP);
        // printf(
        //     "sspMod: Quad: ray is outside the box where the soundspeed is defined\n");
    }

    UpdateSSPSegment(x.y, t.y, ssp->z, ssp->NPts, iSeg.z);
    UpdateSSPSegment(x.x, t.x, ssp->Seg.r, ssp->Nr, iSeg.r);
    LinInterpDensity(x.y, ssp, iSeg, o.rho);
    if(iSeg.z >= ssp->Nz - 1 || iSeg.r >= ssp->Nr - 1) {
        RunError(errState, BHC_ERR_QUAD_ISEG);
        // printf(
        //     "iSeg error in Quad: z %d/%d r %d/%d\n", iSeg.z, ssp->Nz, iSeg.r,
        //     ssp->Nr);
        iSeg.z = iSeg.r = 0;
    }

    // for this depth, x.y, get the sound speed at both ends of the segment
    int32_t Nr = ssp->Nr;
    cz1        = ssp->czMat[iSeg.z * Nr + iSeg.r];
    cz2        = ssp->czMat[iSeg.z * Nr + iSeg.r + 1];

    s2      = x.y - ssp->z[iSeg.z];
    delta_z = ssp->z[iSeg.z + 1] - ssp->z[iSeg.z];

    c1 = ssp->cMat[iSeg.z * Nr + iSeg.r] + s2 * cz1;
    c2 = ssp->cMat[iSeg.z * Nr + iSeg.r + 1] + s2 * cz2;

    // s1 = proportional distance of x.x in range
    delta_r = ssp->Seg.r[iSeg.r + 1] - ssp->Seg.r[iSeg.r];
    s1      = (x.x - ssp->Seg.r[iSeg.r]) / delta_r;
    s1      = bhc::min(s1, RL(1.0)); // force piecewise constant extrapolation for points
                                     // outside the box
    s1 = bhc::max(s1, RL(0.0));      // "

    real c = (RL(1.0) - s1) * c1 + s1 * c2;

    // interpolate the attenuation !!!! This will use the wrong segment if the ssp in the
    // envil is sampled at different depths
    s2 = s2 / delta_z; // convert to a proportional depth in the layer
    // volume attenuation is taken from the single c(z) profile
    real cimag = ((RL(1.0) - s2) * ssp->c[iSeg.z] + s2 * ssp->c[iSeg.z + 1]).imag();

    o.ccpx = cpx(c, cimag);

    cz = (RL(1.0) - s1) * cz1 + s1 * cz2;

    cr    = (c2 - c1) / delta_r;
    o.crz = (cz2 - cz1) / delta_r;

    o.gradc = vec2(cr, cz);
    o.crr   = RL(0.0);
    o.czz   = RL(0.0);
}

/**
 * Trilinear hexahedral interpolation of SSP data in 3D
 * assumes a rectilinear case (not the most general hexahedral)
 */
HOST_DEVICE inline void Hexahedral(SSP_3D_FN_ARGS)
{
    if(x.x < ssp->Seg.x[0] || x.x > ssp->Seg.x[ssp->Nx - 1] || x.y < ssp->Seg.y[0]
       || x.y > ssp->Seg.y[ssp->Ny - 1]) {
        RunError(errState, BHC_ERR_OUTSIDE_SSP);
        // printf(
        //     "sspMod: Hexahedral: ray is outside the box where the ocean soundspeed is "
        //     "defined\nx = (x, y, z) = %g, %g, %g\n",
        //     x.x, x.y, x.z);
    }

    UpdateSSPSegment(x.x, t.x, ssp->Seg.x, ssp->Nx, iSeg.x);
    UpdateSSPSegment(x.y, t.y, ssp->Seg.y, ssp->Ny, iSeg.y);
    UpdateSSPSegment(x.z, t.z, ssp->Seg.z, ssp->Nz, iSeg.z);

    // cz at the corners of the current rectangle
    real cz11 = ssp->czMat[((iSeg.x) * ssp->Ny + iSeg.y) * (ssp->Nz - 1) + iSeg.z];
    real cz12 = ssp->czMat[((iSeg.x + 1) * ssp->Ny + iSeg.y) * (ssp->Nz - 1) + iSeg.z];
    real cz21 = ssp->czMat[((iSeg.x) * ssp->Ny + iSeg.y + 1) * (ssp->Nz - 1) + iSeg.z];
    real cz22
        = ssp->czMat[((iSeg.x + 1) * ssp->Ny + iSeg.y + 1) * (ssp->Nz - 1) + iSeg.z];

    // for this depth, x.z get the sound speed at the corners of the current rectangle
    real s3  = x.z - ssp->Seg.z[iSeg.z];
    real c11 = ssp->cMat[((iSeg.x) * ssp->Ny + iSeg.y) * ssp->Nz + iSeg.z] + s3 * cz11;
    real c21 = ssp->cMat[((iSeg.x + 1) * ssp->Ny + iSeg.y) * ssp->Nz + iSeg.z]
        + s3 * cz12;
    real c12 = ssp->cMat[((iSeg.x) * ssp->Ny + iSeg.y + 1) * ssp->Nz + iSeg.z]
        + s3 * cz21;
    real c22 = ssp->cMat[((iSeg.x + 1) * ssp->Ny + iSeg.y + 1) * ssp->Nz + iSeg.z]
        + s3 * cz22;

    // s1 = proportional distance of x.x in x
    real s1 = (x.x - ssp->Seg.x[iSeg.x]) / (ssp->Seg.x[iSeg.x + 1] - ssp->Seg.x[iSeg.x]);
    // force piecewise constant extrapolation for points outside the box
    s1 = bhc::max(bhc::min(s1, RL(1.0)), RL(0.0));

    // s2 = proportional distance of x.y in y
    real s2 = (x.y - ssp->Seg.y[iSeg.y]) / (ssp->Seg.y[iSeg.y + 1] - ssp->Seg.y[iSeg.y]);
    // force piecewise constant extrapolation for points outside the box
    s2 = bhc::max(bhc::min(s2, RL(1.0)), RL(0.0));

    // interpolate the soundspeed in the x direction, at the two endpoints in y (top and
    // bottom sides of rectangle)
    real c1 = c11 + s1 * (c21 - c11);
    real c2 = c12 + s1 * (c22 - c12);
    // c = (RL(1.0) - s2) * c1 + s2 * c2; // interpolation in y
    real cy = (c2 - c1) / (ssp->Seg.y[iSeg.y + 1] - ssp->Seg.y[iSeg.y]);

    // interpolate the soundspeed in the y direction, at the two endpoints in x (left and
    // right sides of rectangle)
    c1 = c11 + s2 * (c12 - c11);
    c2 = c21 + s2 * (c22 - c21);

    real c = c1 + s1 * (c2 - c1); // interpolation in x

    // interpolate the attenuation
    //!! This will use the wrong segment if the ssp in the envfil
    //   is sampled at different depths [LP: this uses ssp->z rather than
    //   ssp->Seg.z]
    s3 /= ssp->z[iSeg.z + 1] - ssp->z[iSeg.z]; // convert s3 to a proportional distance in
                                               // the layer
    // volume attenuation is taken from the single c(z) profile
    real cimag = ((RL(1.0) - s3) * ssp->c[iSeg.z] + s3 * ssp->c[iSeg.z + 1]).imag();
    o.ccpx     = cpx(c, cimag);

    real cx = (c2 - c1) / (ssp->Seg.x[iSeg.x + 1] - ssp->Seg.x[iSeg.x]);

    // same thing on cz
    real cz1 = cz11 + s2 * (cz21 - cz11);
    real cz2 = cz12 + s2 * (cz22 - cz12);
    real cz  = cz1 + s1 * (cz2 - cz1); // interpolation in z

    // o.gradc = vec3(cx, cy, cz);
    o.gradc.x = cx;
    o.gradc.y = cy;
    o.gradc.z = cz;

    o.cxx = RL(0.0);
    o.cyy = RL(0.0);
    o.czz = RL(0.0);
    o.cxy = RL(0.0);
    o.cxz = RL(0.0);
    o.cyz = RL(0.0);

    // linear interpolation for density
    LinInterpDensity(x.z, ssp, iSeg, o.rho);
}

template<bool O3D> HOST_DEVICE inline void Analytic(SSP_TEMPL_FN_ARGS)
{
    iSeg.z  = 0;
    real c0 = FL(1500.0);
    o.rho   = FL(1.0);

    const float unk1 = FL(1300.0);
    const float unk2 = FL(0.00737);

    if constexpr(O3D) {
        real w;
        // if(x.z < 5000.0){
        const float unk3 = FL(100000.0);
        const float unk4 = FL(0.003);
        real epsilon     = unk2 + x.y / unk3 * unk4;
        real epsilon_y   = unk4 / unk3;

        w       = FL(2.0) * (x.z - unk1) / unk1;
        real wz = FL(2.0) / unk1;

        real emw  = STD::exp(-w);
        o.ccpx    = cpx(c0 * (FL(1.0) + epsilon * (w - FL(1.0) + emw)), FL(0.0));
        o.gradc.y = c0 * epsilon_y * (w - FL(1.0) + emw);
        o.gradc.z = c0 * epsilon * (FL(1.0) - emw) * wz;
        o.czz     = c0 * epsilon * emw * SQ(wz);
        o.cyz     = c0 * epsilon_y * (FL(1.0) - emw) * wz;
        // else{ // HOMOGENEOUS HALF-SPACE
        //     w      = FL(2.0) * (FL(5000.0) - unk1) / unk1;
        //     o.ccpx = cpx(c0 * (FL(1.0) + unk2 * (w - FL(1.0) + STD::exp(-w))), FL(0.0);
        //     o.gradc.y = FL(0.0);
        //     o.gradc.z = FL(0.0);
        //     o.czz     = FL(0.0);
        //     o.cyz     = FL(0.0);
        // }

        o.gradc.x = FL(0.0);
        o.cxx     = FL(0.0);
        o.cyy     = FL(0.0);
        o.cxz     = FL(0.0);
        o.cxy     = FL(0.0);
    } else {
        real cr, cz, DxtDz, xt;

        // homogeneous halfspace was removed since BELLHOP needs to get gradc just a
        // little below the boundaries, on ray reflection

        // if(x.y < 5000.0){
        xt        = FL(2.0) * (x.y - unk1) / unk1;
        real emxt = STD::exp(-xt);
        DxtDz     = FL(2.0) / unk1;
        o.ccpx    = cpx(c0 * (FL(1.0) + unk2 * (xt - FL(1.0) + emxt)), FL(0.0));
        cz        = c0 * unk2 * (FL(1.0) - emxt) * DxtDz;
        o.czz     = c0 * unk2 * emxt * SQ(DxtDz);
        //}else{
        // Homogeneous half-space
        // xt     = FL(2.0) * (FL(5000.0) - unk1) / unk1;
        // o.ccpx = cpx(c0 * (FL(1.0) + unk2 * (xt - FL(1.0) + emxt)), FL(0.0));
        // cz     = FL(0.0);
        // o.czz  = FL(0.0);
        //}

        cr      = FL(0.0);
        o.gradc = vec2(cr, cz);
        o.crz   = FL(0.0);
        o.crr   = FL(0.0);
    }
}

template<typename CFG, bool O3D, bool R3D> HOST_DEVICE inline void EvaluateSSP(
    const VEC23<R3D> &x, const VEC23<R3D> &t, SSPOutputs<R3D> &o,
    const Origin<O3D, R3D> &org, const SSPStructure *ssp, SSPSegState &iSeg,
    ErrState *errState)
{
    VEC23<O3D> x_proc = RayToOceanX(x, org);
    VEC23<O3D> t_proc = RayToOceanT(t, org);
    SSPOutputs<O3D> o_proc;
    if constexpr(CFG::ssp::Is1D()) {
        vec2 x_rz, t_rz;
        SSPOutputs<false> o_rz;
        if constexpr(O3D) {
            x_rz = vec2(RL(0.0), x_proc.z);
            t_rz = vec2(RL(0.0), t_proc.z);
        } else {
            x_rz = x_proc;
            t_rz = t_proc;
        }
        if constexpr(CFG::ssp::IsN2Linear()) { // N2-linear profile option
            n2Linear(x_rz, t_rz, o_rz, ssp, iSeg, errState);
        } else if constexpr(CFG::ssp::IsCLinear()) { // C-linear profile option
            cLinear(x_rz, t_rz, o_rz, ssp, iSeg, errState);
        } else if constexpr(CFG::ssp::IsCCubic()) { // Cubic spline profile option
            cCubic(x_rz, t_rz, o_rz, ssp, iSeg, errState);
        } else if constexpr(CFG::ssp::IsCPCHIP()) { // monotone PCHIP ACS profile option
            cPCHIP(x_rz, t_rz, o_rz, ssp, iSeg, errState);
        } else {
            static_assert(!sizeof(CFG), "Invalid template in EvaluateSSP");
        }
        if constexpr(O3D) {
            o_proc.gradc = vec3(RL(0.0), RL(0.0), o_rz.gradc.y);
            o_proc.cxx = o_proc.cyy = o_proc.cxy = o_proc.cxz = o_proc.cyz = RL(0.0);
            o_proc.czz                                                     = o_rz.czz;
            o_proc.ccpx                                                    = o_rz.ccpx;
            o_proc.rho                                                     = o_rz.rho;
        } else {
            o_proc = o_rz;
        }
    } else if constexpr(CFG::ssp::Is2D()) {
        if constexpr(O3D) {
            // Should already have been checked in SSP::Validate
            RunError(errState, BHC_ERR_TEMPLATE);
            o_proc.ccpx  = cpx(NAN, NAN);
            o_proc.gradc = vec3(NAN, NAN, NAN);
            o_proc.rho = o_proc.czz = NAN;
            o_proc.cxx = o_proc.cyy = o_proc.cxy = o_proc.cxz = o_proc.cyz = NAN;
        } else {
            if constexpr(CFG::ssp::IsQuad()) {
                Quad(x_proc, t_proc, o_proc, ssp, iSeg, errState);
            }
        }
    } else if constexpr(CFG::ssp::Is3D()) {
        if constexpr(!O3D) {
            RunError(errState, BHC_ERR_TEMPLATE);
            o_proc.ccpx  = cpx(NAN, NAN);
            o_proc.gradc = vec2(NAN, NAN);
            o_proc.rho = o_proc.czz = NAN;
            o_proc.crr = o_proc.crz = NAN;
        } else {
            if constexpr(CFG::ssp::IsHexahedral()) {
                Hexahedral(x_proc, t_proc, o_proc, ssp, iSeg, errState);
            }
        }
    } else if constexpr(CFG::ssp::IsAnyD()) {
        if constexpr(CFG::ssp::IsAnalytic()) { // Analytic profile option
            Analytic<O3D>(x_proc, t_proc, o_proc, ssp, iSeg, errState);
        } else {
            static_assert(!sizeof(CFG), "Invalid template in EvaluateSSP");
        }
    } else {
        static_assert(!sizeof(CFG), "Invalid template in EvaluateSSP");
    }
    if constexpr(O3D && !R3D) {
        o.gradc.x
            = glm::dot(org.tradial, vec2(o_proc.gradc.x, o_proc.gradc.y)); // r derivative
        o.gradc.y = o_proc.gradc.z;                                        // z derivative

        o.crz = org.tradial.x * o_proc.cxz + org.tradial.y * o_proc.cyz;
        o.crr = o_proc.cxx * SQ(org.tradial.x)
            + FL(2.0) * o_proc.cxy * org.tradial.x * org.tradial.y
            + o_proc.cyy * SQ(org.tradial.y);

        o.ccpx = o_proc.ccpx;
        o.rho  = o_proc.rho;
        o.czz  = o_proc.czz;
    } else {
        o = o_proc;
    }
}

} // namespace bhc
