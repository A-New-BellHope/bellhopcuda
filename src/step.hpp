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
#include "boundary.hpp"
#include "ssp.hpp"

namespace bhc {

// #define STEP_DEBUGGING 1

#ifdef BHC_USE_FLOATS
#define INFINITESIMAL_STEP_SIZE (RL(1e-3))
#else
#define INFINITESIMAL_STEP_SIZE (RL(1e-6))
#endif

/**
 * interface crossing in depth
 * LP: 3D only:
 * Step reduction is not done for the top or bottom layer
 * Instead the SSP is extrapolated
 * This prevents problems when the boundaries are outside the domain of the SSP
 */
template<bool O3D> HOST_DEVICE inline void DepthInterfaceCrossing(
    bool stepTo, real &hInt, VEC23<O3D> &x, const VEC23<O3D> &x0, const VEC23<O3D> &urayt,
    const SSPSegState &iSeg0, const SSPStructure *ssp, int32_t &snapDim)
{
    if(!stepTo) hInt = REAL_MAX;
    if(STD::abs(DEP(urayt)) > REAL_EPSILON) {
        if(ssp->z[iSeg0.z] > DEP(x) && (!O3D || iSeg0.z > 0)) {
            hInt = (ssp->z[iSeg0.z] - DEP(x0)) / DEP(urayt);
#ifdef STEP_DEBUGGING
            printf(
                "Shallower bound SSP Z %g > z %g; hInt = %g\n", ssp->z[iSeg0.z], DEP(x),
                hInt);
#endif
            if(stepTo) {
                x       = x0 + hInt * urayt; // X or X,Y
                DEP(x)  = ssp->z[iSeg0.z];
                snapDim = ZDIM<O3D>();
#ifdef STEP_DEBUGGING
                printf("to %20.17f %20.17f\n", x.x, x.y);
#endif
            }
        } else if(ssp->z[iSeg0.z + 1] < DEP(x) && (!O3D || iSeg0.z + 1 < ssp->Nz - 1)) {
            hInt = (ssp->z[iSeg0.z + 1] - DEP(x0)) / DEP(urayt);
#ifdef STEP_DEBUGGING
            printf(
                "Deeper bound SSP Z %g < z %g; hInt = %g\n", ssp->z[iSeg0.z + 1], DEP(x),
                hInt);
#endif
            if(stepTo) {
                x       = x0 + hInt * urayt; // X or X,Y
                DEP(x)  = ssp->z[iSeg0.z + 1];
                snapDim = ZDIM<O3D>();
#ifdef STEP_DEBUGGING
                printf("to %20.17f %20.17f\n", x.x, x.y);
#endif
            }
        }
    }
}

/**
 * ray mask using a box centered at (0, 0) [2D] / (source x, source y, 0) [3D]
 */
template<bool O3D, int DIM> HOST_DEVICE inline void BeamBoxCrossing(
    bool stepTo, real &hBox, VEC23<O3D> &x, const VEC23<O3D> &x0, const VEC23<O3D> &urayt,
    const BeamStructure<O3D> *Beam, const VEC23<O3D> &xs, int32_t &snapDim)
{
    if(!stepTo) hBox = REAL_MAX;
    if(IsOutsideBeamBoxDim<O3D, DIM>(x, Beam, xs)) {
        real d0 = x0[DIM] - BeamBoxCenter<O3D>(xs)[DIM];
        hBox    = (Beam->Box[DIM] - STD::abs(d0)) / STD::abs(urayt[DIM]);
#ifdef STEP_DEBUGGING
        printf(
            "Beam box crossing %c %g, hBox = %g\n",
            (O3D
                 ? (DIM == 0       ? 'X'
                        : DIM == 1 ? 'Y'
                                   : 'Z')
                 : (DIM == 0 ? 'R' : 'Z')),
            Beam->Box[DIM], hBox);
#endif
        if(stepTo) {
            x       = x0 + hBox * urayt;
            x[DIM]  = BeamBoxCenter<O3D>(xs)[DIM] + STD::copysign(Beam->Box[DIM], d0);
            snapDim = DIM;
        }
    }
}

/**
 * LP: h = hTop or hBot
 */
template<bool O3D> HOST_DEVICE inline void TopBotCrossing(
    bool stepTo, real &h, const BdryStateTopBot<O3D> &bd, VEC23<O3D> &x,
    const VEC23<O3D> &x0, const VEC23<O3D> &urayt, bool &refl, int32_t &snapDim)
{
    if(!stepTo) h = REAL_MAX;
    VEC23<O3D> d, d0;
    d  = x - bd.x;  // vector from top / bottom to ray
    d0 = x0 - bd.x; // vector from top / bottom node to ray origin
    // Originally, this value had to be > a small positive number, meaning the
    // new step really had to be outside the boundary, not just to the boundary.
    // Also, this is not missing a normalization factor, bd.n is normalized so
    // this is actually the distance above the top / below the bottom in meters.
    real w = glm::dot(bd.n, d);
    if(stepTo ? (w > -INFINITESIMAL_STEP_SIZE) : (w >= RL(0.0))) {
        h = -glm::dot(d0, bd.n) / glm::dot(urayt, bd.n);
#ifdef STEP_DEBUGGING
        printf("Top/bot crossing h %g\n", h);
#endif
        if(stepTo) {
            x = x0 + h * urayt;
            // Snap to exact top / bot depth value if it's flat
            if(STD::abs(bd.n.x) < REAL_EPSILON
               && (!O3D || STD::abs(bd.n.y) < REAL_EPSILON)) {
                DEP(x) = DEP(bd.x);
            }
            snapDim = ZDIM<O3D>(); // Even if not flat, exactness in Z most important
#ifdef STEP_DEBUGGING
            printf("to %20.17f %20.17f\n", x.x, x.y);
#endif
        }
        refl = true;
    } else {
        refl = false;
    }
}

/**
 * top, bottom, or ocean segment crossing in range / x / y
 * LP: w = range, x, or y
 */
template<bool O3D> HOST_DEVICE inline void TopBotSegCrossing(
    bool stepTo, bool isY, real &hSeg, const BdryLimits &TopSeg, const BdryLimits &BotSeg,
    const real *seg_w, int32_t iSeg, VEC23<O3D> &x, const VEC23<O3D> &x0,
    const VEC23<O3D> &urayt, char qh, const SSPStructure *ssp, bool &topRefl,
    bool &botRefl, int32_t &snapDim)
{
    BdryLimits segLim;
    segLim.min = bhc::max(TopSeg.min, BotSeg.min);
    segLim.max = bhc::min(TopSeg.max, BotSeg.max);

    if(ssp->Type == qh) {
        segLim.min = bhc::max(segLim.min, seg_w[iSeg]);
        segLim.max = bhc::min(segLim.max, seg_w[iSeg + 1]);
    }

    real &x_w      = isY ? x.y : x.x;
    real x0_w      = isY ? x0.y : x0.x;
    real urayt_w   = isY ? urayt.y : urayt.x;
    int32_t snap_w = isY ? 1 : 0;
#ifdef STEP_DEBUGGING
    const char *wlbl = O3D ? (isY ? "Y" : "X") : "R";
#endif
    if(!stepTo) hSeg = REAL_MAX;
    if(STD::abs(urayt_w) > REAL_EPSILON) {
        if(x_w < segLim.min) {
            hSeg = -(x0_w - segLim.min) / urayt_w;
            if(stepTo) {
                x       = x0 + hSeg * urayt;
                x_w     = segLim.min;
                snapDim = snap_w;
                topRefl = botRefl = false;
            }
#ifdef STEP_DEBUGGING
            printf(
                "Min bound SSP %s %g > %s %g; hSeg = %g\n", wlbl, segLim.min, wlbl, x_w,
                hSeg);
#endif
        } else if(x_w > segLim.max) {
            hSeg = -(x0_w - segLim.max) / urayt_w;
            if(stepTo) {
                x       = x0 + hSeg * urayt;
                x_w     = segLim.max;
                snapDim = snap_w;
                topRefl = botRefl = false;
            }
#ifdef STEP_DEBUGGING
            printf(
                "Max bound SSP %s %g < %s %g; hSeg = %g\n", wlbl, segLim.max, wlbl, x_w,
                hSeg);
#endif
        }
    }
}

/**
 * triangle crossing within a top / bottom segment
 */
HOST_DEVICE inline void TriDiagCrossing(
    bool stepTo, real &h, BdryStateTopBot<true> &bd, vec3 &x, const vec3 &x0,
    const vec3 &urayt, bool &topRefl, bool &botRefl, int32_t &snapDim, ErrState *errState)
{
    if(!stepTo) h = REAL_MAX;
    if(bd.td.onEdge) {
        // Was previously on the edge--completely ignore the boundary now
        return;
    }
    vec3 d  = x - bd.xmid;  // vector from top / bottom center to ray end
    vec3 d0 = x0 - bd.xmid; // vector from top / bottom center to ray origin
    vec3 tri_n
        = vec3(-(bd.lSeg.y.max - bd.lSeg.y.min), bd.lSeg.x.max - bd.lSeg.x.min, RL(0.0));
    tri_n /= glm::length(tri_n);
    real over_diag_amount = glm::dot(d, tri_n);
    bool newSide          = over_diag_amount >= RL(0.0);
    if(newSide == bd.td.side) {
        // Didn't cross the edge at all
        return;
    }
    bool newOnEdge = STD::abs(over_diag_amount) < TRIDIAG_THRESH;
    if(newOnEdge) {
        // Naturally stepped right to the edge
#ifdef STEP_DEBUGGING
        printf(
            "Naturally stepped to tri diag edge, over_diag_amount %g\n",
            over_diag_amount);
#endif
        return;
    }
    real hnew = -glm::dot(d0, tri_n) / glm::dot(urayt, tri_n);
    if(hnew < RL(0.0)) {
        RunWarning(errState, BHC_WARN_TRIDIAG_H_NEGATIVE);
        h = RL(0.0);
    } else {
        if(stepTo && hnew >= h) { RunWarning(errState, BHC_WARN_TRIDIAG_H_GROWING); }
        h = hnew;
    }
#ifdef STEP_DEBUGGING
    printf(
        "Tri diag crossing h = %g, dot(n, d0) = %g, dot(n, d) = %g\n", h,
        glm::dot(tri_n, d0), over_diag_amount);
#endif
    if(stepTo) {
        x                   = x0 + h * urayt;
        bd.td.justSteppedTo = true;
        bd.td.outgoingSide  = newSide;
        snapDim             = -2; // Snap to X or Y unspecified
        topRefl = botRefl = false;
    }
}

/**
 * calculate a reduced step size, h, that lands on any points where the environment
 * changes
 *
 * x0, urayt: ray coordinate and tangent
 * iSeg0: SSP layer the ray is in
 * Topx, Topn, Botx, Botn: Top, bottom coordinate and normal
 * h: reduced step size
 */
template<bool O3D> HOST_DEVICE inline void ReduceStep(
    const VEC23<O3D> &x0, const VEC23<O3D> &urayt, const SSPSegState &iSeg0,
    BdryState<O3D> &bds, const BeamStructure<O3D> *Beam, const VEC23<O3D> &xs,
    const SSPStructure *ssp, ErrState *errState, real &h, int32_t &iSmallStepCtr)
{
    VEC23<O3D> x;
    real hInt, hBoxxr, hBoxyz, hBoxz_, hTop, hBot, hxSeg, hySeg, hTopDiag, hBotDiag;
    bool dummy;
    int32_t dummy2;

#ifdef STEP_DEBUGGING
    printf("ReduceStep%s\n", O3D ? "3D" : "2D");
#endif

    // Detect SSP interface or boundary crossing and reduce step, if necessary, to land on
    // that crossing. Keep in mind possibility that user put source right on an interface
    // and that multiple events can occur (crossing interface, top, and bottom in a single
    // step).

    x = x0 + h * urayt; // make a trial step

    DepthInterfaceCrossing<O3D>(false, hInt, x, x0, urayt, iSeg0, ssp, dummy2);
    BeamBoxCrossing<O3D, 0>(false, hBoxxr, x, x0, urayt, Beam, xs, dummy2);
    BeamBoxCrossing<O3D, 1>(false, hBoxyz, x, x0, urayt, Beam, xs, dummy2);
    if constexpr(O3D) {
        BeamBoxCrossing<O3D, 2>(false, hBoxz_, x, x0, urayt, Beam, xs, dummy2);
    } else {
        hBoxz_ = REAL_MAX;
    }
    TopBotCrossing<O3D>(false, hTop, bds.top, x, x0, urayt, dummy, dummy2);
    TopBotCrossing<O3D>(false, hBot, bds.bot, x, x0, urayt, dummy, dummy2);

    if constexpr(O3D) {
        TopBotSegCrossing<O3D>(
            false, false, hxSeg, bds.top.lSeg.x, bds.bot.lSeg.x, ssp->Seg.x, iSeg0.x, x,
            x0, urayt, 'H', ssp, dummy, dummy, dummy2);
        TopBotSegCrossing<O3D>(
            false, true, hySeg, bds.top.lSeg.y, bds.bot.lSeg.y, ssp->Seg.y, iSeg0.y, x,
            x0, urayt, 'H', ssp, dummy, dummy, dummy2);
        TriDiagCrossing(
            false, hTopDiag, bds.top, x, x0, urayt, dummy, dummy, dummy2, errState);
        TriDiagCrossing(
            false, hBotDiag, bds.bot, x, x0, urayt, dummy, dummy, dummy2, errState);
    } else {
        TopBotSegCrossing<O3D>(
            false, false, hxSeg, bds.top.lSeg, bds.bot.lSeg, ssp->Seg.r, iSeg0.r, x, x0,
            urayt, 'Q', ssp, dummy, dummy, dummy2);
        hySeg = hTopDiag = hBotDiag = REAL_MAX;
    }

    // take limit set by shortest distance to a crossing
    h = bhc::min(
        bhc::min(
            bhc::min(h, hInt),
            bhc::min(bhc::min(bhc::min(hBoxxr, hBoxyz), hBoxz_), bhc::min(hTop, hBot))),
        bhc::min(bhc::min(hxSeg, hySeg), bhc::min(hTopDiag, hBotDiag)));

    if(h < RL(-1e-4)) {
        RunWarning(errState, BHC_WARN_STEP_NEGATIVE_H);
        // printf("ReduceStep: negative h %f\n", h);
    }
    if(h < INFINITESIMAL_STEP_SIZE * Beam->deltas) { // is it taking an infinitesimal
                                                     // step?
        h = INFINITESIMAL_STEP_SIZE * Beam->deltas;  // make sure we make some motion
        ++iSmallStepCtr; // keep a count of the number of sequential small steps
#ifdef STEP_DEBUGGING
        printf("Small step forced to %g\n", h);
#endif
    } else {
        iSmallStepCtr = 0;
    }
}

/**
 * snapDim: See OceanToRayX.
 */
template<bool O3D> HOST_DEVICE inline void StepToBdry(
    const VEC23<O3D> &x0, VEC23<O3D> &x2, const VEC23<O3D> &urayt, real &h, bool &topRefl,
    bool &botRefl, int32_t &snapDim, const SSPSegState &iSeg0, BdryState<O3D> &bds,
    const BeamStructure<O3D> *Beam, const VEC23<O3D> &xs, const SSPStructure *ssp,
    ErrState *errState)
{
#ifdef STEP_DEBUGGING
    printf("StepToBdry\n");
#endif
    // Original step due to maximum step size
    h       = Beam->deltas;
    x2      = x0 + h * urayt;
    snapDim = -1;

    DepthInterfaceCrossing<O3D>(true, h, x2, x0, urayt, iSeg0, ssp, snapDim);
    BeamBoxCrossing<O3D, 0>(true, h, x2, x0, urayt, Beam, xs, snapDim);
    BeamBoxCrossing<O3D, 1>(true, h, x2, x0, urayt, Beam, xs, snapDim);
    if constexpr(O3D) {
        BeamBoxCrossing<O3D, 2>(true, h, x2, x0, urayt, Beam, xs, snapDim);
    }
    TopBotCrossing<O3D>(true, h, bds.top, x2, x0, urayt, topRefl, snapDim);
    TopBotCrossing<O3D>(true, h, bds.bot, x2, x0, urayt, botRefl, snapDim);
    if(botRefl) topRefl = false;

    if constexpr(O3D) {
        TopBotSegCrossing<O3D>(
            true, false, h, bds.top.lSeg.x, bds.bot.lSeg.x, ssp->Seg.x, iSeg0.x, x2, x0,
            urayt, 'H', ssp, topRefl, botRefl, snapDim);
        TopBotSegCrossing<O3D>(
            true, true, h, bds.top.lSeg.y, bds.bot.lSeg.y, ssp->Seg.y, iSeg0.y, x2, x0,
            urayt, 'H', ssp, topRefl, botRefl, snapDim);
        TriDiagCrossing(
            true, h, bds.top, x2, x0, urayt, topRefl, botRefl, snapDim, errState);
        TriDiagCrossing(
            true, h, bds.bot, x2, x0, urayt, topRefl, botRefl, snapDim, errState);
    } else {
        TopBotSegCrossing<O3D>(
            true, false, h, bds.top.lSeg, bds.bot.lSeg, ssp->Seg.r, iSeg0.r, x2, x0,
            urayt, 'Q', ssp, topRefl, botRefl, snapDim);
    }

    // is it taking an infinitesimal step?
    if(h < INFINITESIMAL_STEP_SIZE * Beam->deltas) {
        h       = INFINITESIMAL_STEP_SIZE * Beam->deltas; // make sure we make some motion
        x2      = x0 + h * urayt;
        snapDim = -1;
#ifdef STEP_DEBUGGING
        printf("StepToBdry small step forced h %g to (%g,%g)\n", h, x2.x, x2.y);
#endif
        // Recheck reflection conditions
        VEC23<O3D> d;
        d = x2 - bds.top.x; // vector from top to ray
        if(glm::dot(bds.top.n, d) > REAL_EPSILON) {
            topRefl = true;
            if constexpr(O3D) {
                bds.top.td.justSteppedTo = false;
                bds.bot.td.justSteppedTo = false;
            }
            snapDim = ZDIM<O3D>();
        } else {
            topRefl = false;
        }
        d = x2 - bds.bot.x; // vector from bottom to ray
        if(glm::dot(bds.bot.n, d) > REAL_EPSILON) {
            botRefl = true;
            topRefl = false;
            if constexpr(O3D) {
                bds.top.td.justSteppedTo = false;
                bds.bot.td.justSteppedTo = false;
            }
            snapDim = ZDIM<O3D>();
        } else {
            botRefl = false;
        }
    } else {
#ifdef STEP_DEBUGGING
        printf("StepToBdry normal h %20.17g to (%20.17g,%20.17g)\n", h, x2.x, x2.y);
#endif
    }
}

template<bool R3D> HOST_DEVICE inline void Get_c_partials(
    const rayPt<R3D> &ray, const SSPOutputs<R3D> &o, StepPartials<R3D> &part)
{
    if constexpr(R3D) {
        vec3 e1, e2;
        RayNormal(ray.t, ray.phi, o.ccpx.real(), e1, e2);

        part.cnn = o.cxx * SQ(e1.x) + o.cyy * SQ(e1.y) + o.czz * SQ(e1.z)
            + FL(2.0) * o.cxy * e1.x * e1.y + FL(2.0) * o.cxz * e1.x * e1.z
            + FL(2.0) * o.cyz * e1.y * e1.z;

        part.cmn = o.cxx * e1.x * e2.x + o.cyy * e1.y * e2.y + o.czz * e1.z * e2.z
            + o.cxy * (e1.x * e2.y + e2.x * e1.y) + o.cxz * (e1.x * e2.z + e2.x * e1.z)
            + o.cyz * (e1.y * e2.z + e2.y * e1.z);

        part.cmm = o.cxx * SQ(e2.x) + o.cyy * SQ(e2.y) + o.czz * SQ(e2.z)
            + FL(2.0) * o.cxy * e2.x * e2.y + FL(2.0) * o.cxz * e2.x * e2.z
            + FL(2.0) * o.cyz * e2.y * e2.z;
    } else {
        part.cnn_csq = o.crr * SQ(ray.t.y) - FL(2.0) * o.crz * ray.t.x * ray.t.y
            + o.czz * SQ(ray.t.x);
    }
}

template<bool R3D> HOST_DEVICE inline rayPtExtras<R3D> ComputeDeltaPQ(
    const rayPt<R3D> &ray, const SSPOutputs<R3D> &o, const StepPartials<R3D> &part)
{
    rayPtExtras<R3D> pq;
    if constexpr(R3D) {
        pq.phi = (FL(1.0) / o.ccpx.real()) * ray.t.z
            * (ray.t.y * o.gradc.x - ray.t.x * o.gradc.y) / (SQ(ray.t.x) + SQ(ray.t.y));

        mat2x2 c_mat(part.cnn, part.cmn, part.cmn, part.cmm);
        c_mat *= -RL(1.0) / SQ(o.ccpx.real());

        // LP: Want to separately multiply _tilde and _hat by c_mat. If they
        // were column vectors, we could just multiply the matrix like normal,
        // but they are row vectors. P^T = C * Q^T --> P = Q * C^T
        pq.p = ray.q * glm::transpose(c_mat);
    } else {
        pq.p = -part.cnn_csq * ray.q;
    }
    pq.q = o.ccpx.real() * ray.p;
    return pq;
}

template<bool R3D> HOST_DEVICE inline void UpdateRayPQ(
    rayPt<R3D> &ray1, const rayPt<R3D> &ray0, real h, const rayPtExtras<R3D> &pq)
{
    if constexpr(R3D) { ray1.phi = ray0.phi + h * pq.phi; }
    ray1.p = ray0.p + h * pq.p;
    ray1.q = ray0.q + h * pq.q;
}

/**
 * LP: Set DMat to zeros to get non-reflect version.
 */
template<bool REFLECTVERSION> HOST_DEVICE inline void CurvatureCorrection3D(
    rayPt<true> &ray, const mat2x2 &DMat, real Tg, real Th, real cn1jump, real cn2jump,
    real csjump, const vec3 &rayn1, const vec3 &rayn2, const vec3 &e1, const vec3 &e2)
{
    real rm = Tg / Th; // this is tan( alpha ) where alpha is the angle of incidence

    // Note that Tg, Th need to be multiplied by c to normalize tangent; hence, csq below
    // added the copysign in r2 to make ati and bty have a symmetric effect on the beam
    // not clear why that's needed

    mat2x2 rmat; // originally R1, R2, R3
    real csq   = SQ(ray.c);
    rmat[0][0] = FL(2.0) / csq * DMat[0][0] / Th
        + rm * (FL(2.0) * cn1jump - rm * csjump) / csq;
    rmat[0][1] = FL(2.0) / ray.c * DMat[1][0] * STD::copysign(RL(1.0), -Th)
        + rm * cn2jump / csq;
    rmat[1][0] = rmat[0][1];
    rmat[1][1] = FL(2.0) * DMat[1][1] * Th;
    // PrintMatrix(rmat, "rmat");

    if constexpr(REFLECTVERSION) {
        rmat[1][0] = -rmat[1][0]; // LP: only first row, second column is negative
        // rmat[1][1] = -rmat[1][1]; // mbp: this one good [LP: The equivalent of this is
        // commented out with this comment]
    } else {
        rmat = -rmat; // LP: all terms are negative
    }

    // z-component of unit tangent is sin( theta ); we want cos( theta )
    // rmat[0][0] *= (FL(1.0) - SQ(ray.c * ray.t.z))

    bool noCurvatureChange = false; // Silence MSVC warning
    if constexpr(REFLECTVERSION) {
        noCurvatureChange = rmat[0][0] == RL(0.0) && rmat[0][1] == RL(0.0)
            && rmat[1][1] == RL(0.0);
    }
    if(noCurvatureChange) {
        // LP: There is no curvature change, but rotating p forward and back
        // can change it slightly due to floating-point imprecision, leading to
        // long-term divergence.
    } else {
        // *** curvature correction ***

        /*
        LP: Arrays in Fortran are stored as A[row][col] where the leftmost index
        (row) is the small increment to adjacent memory. Arrays in C are stored
        as A[row][col] where the rightmost index (col) is the small increment to
        adjacent memory. Arrays in OpenGL and therefore GLM are stored as
        A[col][row] where the rightmost index (row) is the small increment to
        adjacent memory. For bellhopcxx/bellhopcuda, we don't care how the data
        is stored for proper matrices like this, but the indexing has to be
        swapped compared to the Fortran.
        */
        mat2x2 RotMat;
        RotMat[0][0] = glm::dot(rayn1, e1);
        RotMat[1][0] = glm::dot(rayn1, e2);
        RotMat[0][1] = -RotMat[1][0]; // mbp: same as glm::dot(rayn2, e1)
        RotMat[1][1] = glm::dot(rayn2, e2);

        // rotate p-q values in e1, e2 system, onto rayn1, rayn2 system
        // LP: _tilde and _hat are the first and second ROWS of p / q.

        mat2x2 p_in = RotMat * ray.p;
        mat2x2 q_in = RotMat * ray.q;

        // here's the actual curvature change

        mat2x2 p_out = p_in + rmat * q_in;

        // rotate p back to e1, e2 system, q does not change
        // Note RotMat^(-1) = RotMat^T

        ray.p = glm::transpose(RotMat) * p_out;
    }

    if constexpr(REFLECTVERSION) {
        // Logic below fixes a bug when the |dot product| is infinitesimally greater than
        // 1 (then ACos is complex)
        real d;
        d = glm::dot(rayn1, e1);
        // d = rayn1.y * e1.y;
        // d = STD::fma(rayn1.x, e1.x, d);
        // d = STD::fma(rayn1.z, e1.z, d);
        ray.phi += FL(2.0)
            * STD::acos(bhc::max(bhc::min(d, RL(1.0)), RL(-1.0))); // What happens to
                                                                   // torsion?
        // printf("dot %17.12g phi %17.12g\n",
        //     glm::dot(rayn1, e1), ray.phi);
    }
}

HOST_DEVICE inline void CalcTangent_Normals(
    const rayPt<true> &ray, real c, const vec3 &nBdry, vec3 &rayt, vec3 &rayn1,
    vec3 &rayn2, real rayn2sign)
{
    rayt  = c * ray.t;                           // unit tangent to ray
    rayn2 = rayn2sign * glm::cross(rayt, nBdry); // ray tangent x boundary normal gives
                                                 // refl. plane normal
    rayn2 /= glm::length(rayn2);                 // unit normal
    rayn1 = -glm::cross(rayt, rayn2); // ray tangent x refl. plane normal is first ray
                                      // normal
}

/**
 * correct p-q due to jumps in the gradient of the sound speed
 */
template<bool R3D> HOST_DEVICE inline void CurvatureCorrection(
    rayPt<R3D> &ray, const VEC23<R3D> &gradcjump, const SSPSegState &iSeg,
    const SSPSegState &iSeg0)
{
    if constexpr(R3D) {
        vec3 nBdry(RL(0.0), RL(0.0), RL(0.0));
        // what if we cross iSeg.x, iSeg.y, or iSeg.z at the same time?
        if(iSeg.z != iSeg0.z) {
            nBdry.z = -STD::copysign(RL(1.0), ray.t.z); // inward normal to layer
        } else if(iSeg.x != iSeg0.x) {
            nBdry.x = -STD::copysign(RL(1.0), ray.t.x); // inward normal to x-segment
        } else {
            nBdry.y = -STD::copysign(RL(1.0), ray.t.y); // inward normal to y-segment
        }

        real Th = glm::dot(ray.t, nBdry); // component of ray tangent, normal to boundary
        vec3 tBdry = ray.t - Th * nBdry; // tangent, along the boundary, in the reflection
                                         // plane
        tBdry /= glm::length(tBdry);     // unit boundary tangent
        real Tg = glm::dot(ray.t, tBdry); // component of ray tangent, along the boundary

        vec3 rayt, rayn1, rayn2;
        CalcTangent_Normals(ray, ray.c, nBdry, rayt, rayn1, rayn2, RL(1.0));

        vec3 e1, e2;
        RayNormal_unit(rayt, ray.phi, e1, e2);

        // normal and tangential derivatives of the sound speed
        real cn1jump = glm::dot(gradcjump, rayn1);
        real cn2jump = glm::dot(gradcjump, rayn2);
        real csjump  = glm::dot(gradcjump, rayt);

        // printf("cn1 cn2 cs jumps %g %g %g\n", cn1jump, cn2jump, csjump);

        CurvatureCorrection3D<false>(
            ray, mat2x2(RL(0.0)), Tg, Th, cn1jump, cn2jump, csjump, rayn1, rayn2, e1, e2);
    } else {
        // LP: Nx2D only:
        // mbp: this needs modifying like the full 3D version to handle jumps in the x-y
        // direction
        vec2 ray2n = vec2(-ray.t.y, ray.t.x); // ray normal

        real cnjump = glm::dot(gradcjump, ray2n);
        real csjump = glm::dot(gradcjump, ray.t);

        real rm, rn;
        if(iSeg.z != iSeg0.z) {     // crossing in depth
            rm = ray.t.x / ray.t.y; // this is tan( alpha ) where alpha is the angle of
                                    // incidence
        } else {                    // crossing in range
            // LP: This case is excluded for Nx2D by the if condition under
            // which this is called.
            rm = -ray.t.y / ray.t.x; // this is tan( alpha ) where alpha is the angle of
                                     // incidence
        } // LP: The case where it crosses in depth and range simultaneously is not
          // handled.

        rn    = rm * (FL(2.0) * cnjump - rm * csjump) / ray.c;
        ray.p = ray.p - ray.q * rn;
    }
}

/**
 * Does a single step along the ray
 */
template<typename CFG, bool O3D, bool R3D> HOST_DEVICE inline void Step(
    rayPt<R3D> ray0, rayPt<R3D> &ray2, BdryState<O3D> &bds,
    const BeamStructure<O3D> *Beam, const VEC23<O3D> &xs, const Origin<O3D, R3D> &org,
    const SSPStructure *ssp, SSPSegState &iSeg, ErrState *errState,
    int32_t &iSmallStepCtr, bool &topRefl, bool &botRefl)
{
    rayPt<R3D> ray1;
    SSPOutputs<R3D> o0, o1, o2;
    StepPartials<R3D> part0, part1;
    rayPtExtras<R3D> pq0, pq1;
    VEC23<R3D> urayt0, urayt1;
    real csq0, csq1, h, w0, w1, hw0, hw1;

#ifdef STEP_DEBUGGING
    // printf("\nray0 x t p q tau amp (%20.17f,%20.17f) (%20.17f,%20.17f)
    // (%20.17f,%20.17f) (%20.17f,%20.17f) (%20.17f,%20.17f) %20.17f\n",
    //     ray0.x.x, ray0.x.y, ray0.t.x, ray0.t.y, ray0.p.x, ray0.p.y, ray0.q.x, ray0.q.y,
    //     ray0.tau.real(), ray0.tau.imag(), ray0.Amp);
    // printf("iSeg.z iSeg.r %d %d\n", iSeg.z, iSeg.r);
    if constexpr(R3D) {
        printf(
            "\nray0 x t (%20.17f,%20.17f,%20.17f) (%20.17e,%20.17e,%20.17e)\n", ray0.x.x,
            ray0.x.y, ray0.x.z, ray0.t.x, ray0.t.y, ray0.t.z);
        PrintMatrix(ray0.p, "ray0.p");
        PrintMatrix(ray0.q, "ray0.q");
        printf(
            "iSegx iSegy iSegz top.td.side bot.td.side %d %d %d %c %c\n", iSeg.x + 1,
            iSeg.y + 1, iSeg.z + 1, bds.top.td.side ? 'T' : 'F',
            bds.bot.td.side ? 'T' : 'F');
    } else {
        printf(
            "\nray0 x t (%20.17f,%20.17f) (%20.17e,%20.17e)\n", ray0.x.x, ray0.x.y,
            ray0.t.x, ray0.t.y);
        printf("iSegr iSegz %d %d\n", iSeg.r + 1, iSeg.z + 1);
    }
// if(ray0.x.x > RL(10.0)){
//     printf("Enough\n");
//     bail();
// }
#endif

    // The numerical integrator used here is a version of the polygon (a.k.a. midpoint,
    // leapfrog, or Box method), and similar to the Heun (second order Runge-Kutta
    // method). However, it's modified to allow for a dynamic step change, while
    // preserving the second-order accuracy).

    // *** Phase 1 (an Euler step)

    EvaluateSSP<CFG, O3D, R3D>(ray0.x, ray0.t, o0, org, ssp, iSeg, errState);
    // printf("iSeg.z iSeg.r %d %d\n", iSeg.z, iSeg.r);
    Get_c_partials<R3D>(ray0, o0, part0);
    pq0 = ComputeDeltaPQ<R3D>(ray0, o0, part0);

    SSPSegState iSeg0 = iSeg; // make note of current layer

    csq0   = SQ(o0.ccpx.real());
    urayt0 = o0.ccpx.real() * ray0.t; // unit tangent
    h      = Beam->deltas; // initially set the step h, to the basic one, deltas

    // printf("urayt0 (%g,%g)\n", urayt0.x, urayt0.y);

    // reduce h to land on boundary
    VEC23<O3D> x_o = RayToOceanX(ray0.x, org);
    VEC23<O3D> t_o = RayToOceanT(urayt0, org);
    ReduceStep<O3D>(x_o, t_o, iSeg0, bds, Beam, xs, ssp, errState, h, iSmallStepCtr);
    // printf("out h, urayt0 %20.17f (%20.17f, %20.17f)\n", h, urayt0.x, urayt0.y);
    real halfh = FL(0.5) * h; // first step of the modified polygon method is a half step

    ray1.x = ray0.x + halfh * urayt0;
    ray1.t = ray0.t - halfh * o0.gradc / csq0;
    UpdateRayPQ<R3D>(ray1, ray0, halfh, pq0);

    // printf("ray1 x t p q (%20.17f,%20.17f) (%20.17f,%20.17f) (%20.17f,%20.17f)
    // (%20.17f,%20.17f)\n",
    //     ray1.x.x, ray1.x.y, ray1.t.x, ray1.t.y, ray1.p.x, ray1.p.y, ray1.q.x,
    //     ray1.q.y);

    // *** Phase 2

    EvaluateSSP<CFG, O3D, R3D>(ray1.x, ray1.t, o1, org, ssp, iSeg, errState);
    Get_c_partials<R3D>(ray1, o1, part1);
    pq1 = ComputeDeltaPQ<R3D>(ray1, o1, part1);

    // The Munk test case with a horizontally launched ray caused problems.
    // The ray vertexes on an interface and can ping-pong around that interface.
    // Have to be careful in that case about big changes to the stepsize (that invalidate
    // the leap-frog scheme) in phase II. A modified Heun or Box method could also work.

    csq1   = SQ(o1.ccpx.real());
    urayt1 = o1.ccpx.real() * ray1.t; // unit tangent

    // printf("urayt1 (%g,%g)\n", urayt1.x, urayt1.y);

    // reduce h to land on boundary
    t_o = RayToOceanT(urayt1, org);
    ReduceStep<O3D>(x_o, t_o, iSeg0, bds, Beam, xs, ssp, errState, h, iSmallStepCtr);

    // use blend of f' based on proportion of a full step used.
    w1 = h / (RL(2.0) * halfh);
    w0 = RL(1.0) - w1;
    // printf("w1 %20.17f w0 %20.17f\n", w1, w0);
    VEC23<R3D> urayt2 = w0 * urayt0 + w1 * urayt1;
    // Take the blended ray tangent (urayt2) and find the minimum step size (h)
    // to put this on a boundary, and ensure that the resulting position
    // (ray2.x) gets put precisely on the boundary.
    VEC23<O3D> x2_o;
    t_o = RayToOceanT(urayt2, org);
    int32_t snapDim;
    StepToBdry<O3D>(
        x_o, x2_o, t_o, h, topRefl, botRefl, snapDim, iSeg0, bds, Beam, xs, ssp,
        errState);
    ray2.x = OceanToRayX(x2_o, org, urayt2, snapDim, errState);
#ifdef STEP_DEBUGGING
    if constexpr(O3D && !R3D) {
        VEC23<O3D> x2_o_out = RayToOceanX(ray2.x, org);
        printf(
            "OceanToRayX in (%20.17f,%20.17f,%20.17f)\n"
            "           ==> (%20.17f,%20.17f)\n"
            "           ==> (%20.17f,%20.17f,%20.17f)\n",
            x2_o.x, x2_o.y, x2_o.z, ray2.x.x, ray2.x.y, x2_o_out.x, x2_o_out.y,
            x2_o_out.z);
    }
#endif

    // Update other variables with this new h
    hw0      = h * w0;
    hw1      = h * w1;
    ray2.t   = ray0.t - hw0 * o0.gradc / csq0 - hw1 * o1.gradc / csq1;
    ray2.tau = ray0.tau + hw0 / o0.ccpx + hw1 / o1.ccpx;
    UpdateRayPQ<R3D>(ray2, ray0, hw0, pq0);
    UpdateRayPQ<R3D>(ray2, ray2, hw1, pq1); // Not a typo, accumulating into 2

#ifdef STEP_DEBUGGING
    if constexpr(R3D) {
        printf("ray2.t (%20.17f,%20.17f,%20.17f)\n", ray2.t.x, ray2.t.y, ray2.t.z);
    }
#endif

    ray2.Amp       = ray0.Amp;
    ray2.Phase     = ray0.Phase;
    ray2.NumTopBnc = ray0.NumTopBnc;
    ray2.NumBotBnc = ray0.NumBotBnc;

    // If we crossed an interface, apply jump condition

    EvaluateSSP<CFG, O3D, R3D>(ray2.x, ray2.t, o2, org, ssp, iSeg, errState);
    ray2.c = o2.ccpx.real();

    if(iSeg.z != iSeg0.z || (!R3D && !O3D && iSeg.r != iSeg0.r)
       || (R3D && (iSeg.x != iSeg0.x || iSeg.y != iSeg0.y))) {
        VEC23<R3D> gradcjump = o2.gradc - o0.gradc;
        CurvatureCorrection<R3D>(ray2, gradcjump, iSeg, iSeg0);
    }
    // if constexpr(R3D){
    //     PrintMatrix(ray2.p, "ray2.p");
    //     PrintMatrix(ray2.q, "ray2.q");
    // }
}

} // namespace bhc
