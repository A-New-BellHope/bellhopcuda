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
#include "ssp.hpp"
#include "boundary.hpp"

namespace bhc {

//#define STEP_DEBUGGING 1

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
    real &h, VEC23<O3D> &x, const VEC23<O3D> &x0, const VEC23<O3D> &urayt,
    const SSPSegState &iSeg0, const SSPStructure *ssp, bool stepTo)
{
    if(!stepTo) h = REAL_MAX;
    if(STD::abs(DEP(urayt)) > REAL_EPSILON){
        if(      ssp->z[iSeg0.z]     > DEP(x) && (!O3D || iSeg0.z > 0)){
            h = (ssp->z[iSeg0.z]     - DEP(x0)) / DEP(urayt);
            if(stepTo){
                x = x0 + h * urayt; // X or X,Y
                DEP(x) = ssp->z[iSeg0.z];
                #ifdef STEP_DEBUGGING
                printf("to %20.17f %20.17f\n", x.x, x.y);
                #endif
            }
            #ifdef STEP_DEBUGGING
            printf("Shallower bound SSP Z %g > z %g; h = %g\n", ssp->z[iSeg0.z], DEP(x), h);
            #endif
        }else if(ssp->z[iSeg0.z + 1] < DEP(x) && (!O3D || iSeg0.z + 1 < ssp->Nz-1)){
            h = (ssp->z[iSeg0.z + 1] - DEP(x0)) / DEP(urayt);
            #ifdef STEP_DEBUGGING
            printf("Deeper bound SSP Z %g < z %g; h = %g\n", ssp->z[iSeg0.z+1], DEP(x), h);
            #endif
            if(stepTo){
                x = x0 + h * urayt; // X or X,Y
                DEP(x) = ssp->z[iSeg0.z + 1];
                #ifdef STEP_DEBUGGING
                printf("to %20.17f %20.17f\n", x.x, x.y);
                #endif
            }
        }
    }
}

template<bool O3D> HOST_DEVICE inline void TopBotCrossing(
    real &h, const BdryStateTopBot<O3D> &bd, VEC23<O3D> &x,
    const VEC23<O3D> &x0, const VEC23<O3D> &urayt, bool stepTo, bool &refl)
{
    if(!stepTo) h = REAL_MAX;
    VEC23<O3D> d, d0;
    d  = x  - bd.x; // vector from top / bottom to ray
    d0 = x0 - bd.x; // vector from top / bottom node to ray origin
    // Originally, this value had to be > a small positive number, meaning the
    // new step really had to be outside the boundary, not just to the boundary.
    // Also, this is not missing a normalization factor, Topn is normalized so
    // this is actually the distance above the top in meters.
    real w = glm::dot(bd.n, d);
    if(stepTo ? (w > -INFINITESIMAL_STEP_SIZE) : (w >= RL(0.0))){
        h = -glm::dot(d0, bd.n) / glm::dot(urayt, bd.n);
        #ifdef STEP_DEBUGGING
        printf("Top/bot crossing h %g\n", h);
        #endif
        if(stepTo){
            x = x0 + h * urayt;
            // Snap to exact top / bot depth value if it's flat
            if(STD::abs(bd.n.x) < REAL_EPSILON && (!O3D || STD::abs(bd.n.y) < REAL_EPSILON)){
                DEP(x) = DEP(bd.x);
            }
            #ifdef STEP_DEBUGGING
            printf("to %20.17f %20.17f\n", x.x, x.y);
            #endif
        }
        refl = true;
    }else{
        refl = false;
    }
}

/**
 * top or bottom segment crossing in range / x / y
 */
template<bool O3D> HOST_DEVICE inline void TopBotSegCrossing(
    real &h, const BdryLimits &TopSeg, const BdryLimits &BotSeg,
    const real *seg_w, int32_t iSeg, VEC23<O3D> &x,
    const VEC23<O3D> &x0, const VEC23<O3D> &urayt, char qh, 
    const SSPStructure *ssp, bool stepTo, bool &topRefl, bool &botRefl, bool isY)
{
    BdryLimits rSeg;
    rSeg.min = bhc::max(TopSeg.min, BotSeg.min);
    rSeg.max = bhc::min(TopSeg.max, BotSeg.max);
    
    if(ssp->Type == qh){
        rSeg.min = bhc::max(rSeg.min, seg_w[iSeg  ]);
        rSeg.max = bhc::min(rSeg.max, seg_w[iSeg+1]);
    }
    
    real &x_w    = isY ?     x.y :     x.x;
    real x0_w    = isY ?    x0.y :    x0.x;
    real urayt_w = isY ? urayt.y : urayt.x;
    if(!stepTo) h = REAL_MAX;
    if(STD::abs(urayt_w) > REAL_EPSILON){
        if(x_w < rSeg.min){
            h = -(x0_w - rSeg.min) / urayt_w;
            if(stepTo){
                x = x0 + h * urayt;
                x_w = rSeg.min;
                topRefl = botRefl = false;
            }
            #ifdef STEP_DEBUGGING
            printf("Min bound SSP R %g > r %g; h4 = %g\n", rSeg.min, x_w, h);
            #endif
        }else if(x_w > rSeg.max){
            h = -(x0_w - rSeg.max) / urayt_w;
            if(stepTo){
                x = x0 + h * urayt;
                x_w = rSeg.max;
                topRefl = botRefl = false;
            }
            #ifdef STEP_DEBUGGING
            printf("Max bound SSP R %g < r %g; h4 = %g; top.max %g; bot.max %g; ssp r up %g\n",
                rSeg.max, x_w, h, TopSeg.max, BotSeg.max,
                ssp->Type == qh ? seg_w[iSeg+1] : INFINITY);
            #endif
        }
    }
}

/**
 * triangle crossing within a top / bottom segment
 */
HOST_DEVICE inline void TriDiagCrossing(
    real &h, const BdryStateTopBot<true> &bd,
    vec3 &x, const vec3 &x0, const vec3 &urayt,
    bool stepTo, bool &topRefl, bool &botRefl)
{
    vec3 d     = x  - bd.x; // vector from top / bottom node to ray end
    vec3 d0    = x0 - bd.x; // vector from top / bottom node to ray origin
    vec3 tri_n = vec3(-(bd.lSeg.y.max - bd.lSeg.y.min), bd.lSeg.x.max - bd.lSeg.x.min, RL(0.0));
    
    if(!stepTo) h = REAL_MAX;
    if( (glm::dot(tri_n, d0) > RL(0.0) && glm::dot(tri_n, d) <= RL(0.0)) ||
        (glm::dot(tri_n, d0) < RL(0.0) && glm::dot(tri_n, d) >= RL(0.0)) ){
        h = -glm::dot(d0, tri_n) / glm::dot(urayt, tri_n);
        #ifdef STEP_DEBUGGING
        printf("Tri diag crossing h = %g\n", h);
        #endif
        if(stepTo){
            int32_t i;
            for(i=0; i<100; ++i){
                x = x0 + h * urayt;
                // LP: Since this is not an exact floating-point value to step
                // to, make sure we have stepped over the boundary.
                if( (glm::dot(tri_n, d0) > RL(0.0) && glm::dot(tri_n, x - bd.x) <= RL(0.0)) ||
                    (glm::dot(tri_n, d0) < RL(0.0) && glm::dot(tri_n, x - bd.x) >= RL(0.0)) ){
                    break;
                }
                // Slightly increase h if not.
                h *= FL(1.000001);
            }
            if(i == 100){
                printf("Warning, TriDiagCrossing did not converge\n");
            }
            topRefl = botRefl = false;
        }
    }
}

/**
 * calculate a reduced step size, h, that lands on any points where the environment changes
 * 
 * x0, urayt: ray coordinate and tangent
 * iSeg0: SSP layer the ray is in
 * Topx, Topn, Botx, Botn: Top, bottom coordinate and normal
 * h: reduced step size
 */
template<bool O3D> HOST_DEVICE inline void ReduceStep(
    const VEC23<O3D> &x0, const VEC23<O3D> &urayt,
    const SSPSegState &iSeg0, const BdryState<O3D> &bds,
    const BeamStructure *Beam, const SSPStructure *ssp, 
    real &h, int32_t &iSmallStepCtr)
{
    VEC23<O3D> x;
    real h1, h2, h3, h4, h5, h6, h7;
    bool dummy;
    
    // Detect interface or boundary crossing and reduce step, if necessary, to land on that crossing.
    // Keep in mind possibility that user put source right on an interface
    // and that multiple events can occur (crossing interface, top, and bottom in a single step).

    x = x0 + h * urayt; // make a trial step

    DepthInterfaceCrossing<O3D>(h1, x, x0, urayt, iSeg0, ssp, false);
    TopBotCrossing<O3D>(h2, bds.top, x, x0, urayt, false, dummy);
    TopBotCrossing<O3D>(h3, bds.bot, x, x0, urayt, false, dummy);
    
    if constexpr(O3D){
        TopBotSegCrossing<O3D>(h4, bds.top.lSeg.x, bds.bot.lSeg.x, ssp->Seg.x, iSeg0.x,
            x, x0, urayt, 'H', ssp, false, dummy, dummy, false);
        TopBotSegCrossing<O3D>(h5, bds.top.lSeg.y, bds.bot.lSeg.y, ssp->Seg.y, iSeg0.y,
            x, x0, urayt, 'H', ssp, false, dummy, dummy, true);
        TriDiagCrossing(h6, bds.top, x, x0, urayt, false, dummy, dummy);
        TriDiagCrossing(h7, bds.bot, x, x0, urayt, false, dummy, dummy);
    }else{
        TopBotSegCrossing<O3D>(h4, bds.top.lSeg, bds.bot.lSeg, ssp->Seg.r, iSeg0.r,
            x, x0, urayt, 'Q', ssp, false, dummy, dummy, false);
        h5 = h6 = h7 = REAL_MAX;
    }
    
    //printf("ReduceStep h h1 h2 h3 h4 h5 h6 h7 %g %g %g %g %g %g %g %g\n",
    //    h, h1, h2, h3, h4, h5, h6, h7);
    // take limit set by shortest distance to a crossing
    h = bhc::min(h, bhc::min(h1, bhc::min(h2, bhc::min(h3, bhc::min(h4, 
        bhc::min(h5, bhc::min(h6, h7)))))));
    /*
    if(h == h1){
        printf("Step %g due to Z SSP crossing\n", h);
    }else if(h == h2){
        printf("Step %g due to top crossing\n", h);
    }else if(h == h3){
        printf("Step %g due to bottom crossing\n", h);
    }else if(h == h4){
        printf("Step %g due to R/X SSP crossing\n", h);
    }else if(h == h5){
        printf("Step %g due to Y SSP crossing\n", h);
    }else if(h == h6){
        printf("Step %g due to top tri/diag crossing\n", h);
    }else if(h == h7){
        printf("Step %g due to bot tri/diag crossing\n", h);
    }else{
        printf("Step %g (unchanged)\n", h);
    }
    */
    
    if(h < INFINITESIMAL_STEP_SIZE * Beam->deltas){ // is it taking an infinitesimal step?
        h = INFINITESIMAL_STEP_SIZE * Beam->deltas; // make sure we make some motion
        ++iSmallStepCtr; // keep a count of the number of sequential small steps
        // printf("Small step forced to %g\n", h);
    }else{
        iSmallStepCtr = 0;
    }
}

template<bool O3D> HOST_DEVICE inline void StepToBdry(
    const VEC23<O3D> &x0, VEC23<O3D> &x2, const VEC23<O3D> &urayt,
    real &h, bool &topRefl, bool &botRefl,
    const SSPSegState &iSeg0, const BdryState<O3D> &bds,
    const BeamStructure *Beam, const SSPStructure *ssp)
{
    // Original step due to maximum step size
    h = Beam->deltas;
    x2 = x0 + h * urayt;
    
    DepthInterfaceCrossing<O3D>(h, x2, x0, urayt, iSeg0, ssp, true);
    TopBotCrossing<O3D>(h, bds.top, x2, x0, urayt, true, topRefl);
    TopBotCrossing<O3D>(h, bds.bot, x2, x0, urayt, true, botRefl);
    if(botRefl) topRefl = false;
    
    if constexpr(O3D){
        TopBotSegCrossing<O3D>(h, bds.top.lSeg.x, bds.bot.lSeg.x, ssp->Seg.x, iSeg0.x,
            x2, x0, urayt, 'H', ssp, true, topRefl, botRefl, false);
        TopBotSegCrossing<O3D>(h, bds.top.lSeg.y, bds.bot.lSeg.y, ssp->Seg.y, iSeg0.y,
            x2, x0, urayt, 'H', ssp, true, topRefl, botRefl, true);
        TriDiagCrossing(h, bds.top, x2, x0, urayt, true, topRefl, botRefl);
        TriDiagCrossing(h, bds.bot, x2, x0, urayt, true, topRefl, botRefl);
    }else{
        TopBotSegCrossing<O3D>(h, bds.top.lSeg, bds.bot.lSeg, ssp->Seg.r, iSeg0.r,
            x2, x0, urayt, 'Q', ssp, true, topRefl, botRefl, false);
    }
    
    if(h < INFINITESIMAL_STEP_SIZE * Beam->deltas){ // is it taking an infinitesimal step?
        h = INFINITESIMAL_STEP_SIZE * Beam->deltas; // make sure we make some motion
        x2 = x0 + h * urayt;
        #ifdef STEP_DEBUGGING
        printf("StepToBdry small step forced h %g to (%g,%g)\n", h, x2.x, x2.y);
        #endif
        // Recheck reflection conditions
        VEC23<O3D> d;
        d = x2 - bds.top.x; // vector from top to ray
        if(glm::dot(bds.top.n, d) > REAL_EPSILON){
            topRefl = true;
        }else{
            topRefl = false;
        }
        d = x2 - bds.bot.x; // vector from bottom to ray
        if(glm::dot(bds.bot.n, d) > REAL_EPSILON){
            botRefl = true;
            topRefl = false;
        }else{
            botRefl = false;
        }
    }
}

template<bool R3D> HOST_DEVICE inline void Get_c_partials(
    const rayPt<R3D> &ray, const SSPOutputs<R3D> &o,
    StepPartials<R3D> &part)
{
    if constexpr(R3D){
        vec3 e1, e2;
        RayNormal(ray.t, ray.phi, o.ccpx.real(), e1, e2);
        
        part.cnn = o.cxx * SQ(e1.x) + o.cyy * SQ(e1.y) + o.czz * SQ(e1.z)
            + FL(2.0) * o.cxy * e1.x * e1.y
            + FL(2.0) * o.cxz * e1.x * e1.z
            + FL(2.0) * o.cyz * e1.y * e1.z;
        
        part.cmn = o.cxx * e1.x * e2.x + o.cyy * e1.y * e2.y + o.czz * e1.z * e2.z
            + o.cxy * (e1.x * e2.y + e2.x * e1.y)
            + o.cxz * (e1.x * e2.z + e2.x * e1.z)
            + o.cyz * (e1.y * e2.z + e2.y * e1.z);
        
        part.cmm = o.cxx * SQ(e2.x) + o.cyy * SQ(e2.y) + o.czz * SQ(e2.z)
            + FL(2.0) * o.cxy * e2.x * e2.y
            + FL(2.0) * o.cxz * e2.x * e2.z
            + FL(2.0) * o.cyz * e2.y * e2.z;
    }else{
        part.cnn_csq = o.crr * SQ(ray.t.y) 
                    - FL(2.0) * o.crz * ray.t.x * ray.t.y 
                    + o.czz * SQ(ray.t.x);
    }
}

template<bool R3D> HOST_DEVICE inline rayPtExtras<R3D> ComputeDeltaPQ(
    const rayPt<R3D> &ray, const SSPOutputs<R3D> &o,
    const StepPartials<R3D> &part)
{
    rayPtExtras<R3D> pq;
    if constexpr(R3D){
        pq.phi = (FL(1.0) / o.ccpx.real()) * ray.t.z *
            (ray.t.y * o.gradc.x - ray.t.x * o.gradc.y) /
            (SQ(ray.t.x) + SQ(ray.t.y));
        
        mat2x2 c_mat(part.cnn, part.cmn, part.cmn, part.cmm);
        c_mat *= -RL(1.0) / SQ(o.ccpx.real());
        
        pq.p_tilde = c_mat         * ray.q_tilde;
        pq.q_tilde = o.ccpx.real() * ray.p_tilde;
        
        pq.p_hat   = c_mat         * ray.q_hat;
        pq.q_hat   = o.ccpx.real() * ray.p_hat;
        
    }else{
        pq.p = -part.cnn_csq * ray.q;
        pq.q = o.ccpx.real() * ray.p;
    }
    return pq;
}

template<bool R3D> HOST_DEVICE inline void UpdateRayPQ(
    rayPt<R3D> &ray1, const rayPt<R3D> &ray0, real h, const rayPtExtras<R3D> &pq)
{
    if constexpr(R3D){
        ray1.phi     = ray0.phi     + h * pq.phi;
        ray1.p_tilde = ray0.p_tilde + h * pq.p_tilde;
        ray1.q_tilde = ray0.q_tilde + h * pq.q_tilde;
        ray1.p_hat   = ray0.p_hat   + h * pq.p_hat;
        ray1.q_hat   = ray0.q_hat   + h * pq.q_hat;
    }else{
        ray1.p = ray0.p + h * pq.p;
        ray1.q = ray0.q + h * pq.q;
    }
}

/**
 * correct p-q due to jumps in the gradient of the sound speed
 */
template<bool R3D> HOST_DEVICE inline void CurvatureCorrection(
    rayPt<R3D> &ray2, const VEC23<R3D> &gradcjump,
    const SSPSegState &iSeg, const SSPSegState &iSeg0)
{
    if constexpr(R3D){
        vec3 nBdry(RL(0.0), RL(0.0), RL(0.0));
        // what if we cross iSeg.x, iSeg.y, or iSeg.z at the same time?
        if(iSeg.z != iSeg0.z){
            nBdry.z = -STD::copysign(RL(1.0), ray2.t.z); // inward normal to layer
        }else if(iSeg.x != iSeg0.x){
            nBdry.x = -STD::copysign(RL(1.0), ray2.t.x); // inward normal to x-segment
        }else{
            nBdry.y = -STD::copysign(RL(1.0), ray2.t.y); // inward normal to y-segment
        }
        
        real Th    = glm::dot(ray2.t, nBdry); // component of ray tangent, normal to boundary
        vec3 tBdry = ray2.t - Th * nBdry;     // tangent, along the boundary, in the reflection plane
        tBdry     /= glm::length(tBdry);      // unit boundary tangent
        real Tg    = glm::dot(ray2.t, tBdry); // component of ray tangent, along the boundary

        vec3 rayt  = ray2.c * ray2.t;         // unit tangent to ray

        vec3 rayn2 = glm::cross(rayt, nBdry); // ray tangent x boundary normal gives refl. plane normal
        rayn2 /= glm::length(rayn2);          // unit normal
        vec3 rayn1 = -glm::cross(rayt, rayn2);// ray tangent x refl. plane normal is first ray normal
        
        // normal and tangential derivatives of the sound speed
        real cn1jump = glm::dot(gradcjump, rayn1);
        real cn2jump = glm::dot(gradcjump, rayn2);
        real csjump  = glm::dot(gradcjump, rayt);
        
        real rm = Tg / Th; // this is tan( alpha ) where alpha is the angle of incidence
        real r1 = rm * (FL(2.0) * cn1jump - rm * csjump) / SQ(ray2.c);
        real r2 = rm * cn2jump / SQ(ray2.c);
        
        // *** curvature correction ***
        
        vec3 e1, e2;
        RayNormal_unit(rayt, ray2.phi, e1, e2); // Compute ray normals e1 and e2
        
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
        mat2x2 RotMat, pmat, qmat;
        RotMat[0][0] = glm::dot(rayn1, e1);
        RotMat[1][0] = glm::dot(rayn1, e2);
        RotMat[0][1] = -RotMat[1][0]; // glm::dot(rayn2, e1)
        RotMat[1][1] = glm::dot(rayn2, e2);
        
        // rotate p-q values in e1, e2 system, onto rayn1, rayn2 system
        
        glm::row(pmat, 0, ray2.p_tilde);
        glm::row(pmat, 1, ray2.p_hat);
        pmat = RotMat * pmat;
        vec2 p_tilde_in = glm::row(pmat, 0);
        vec2 p_hat_in   = glm::row(pmat, 1);
        
        glm::row(qmat, 0, ray2.q_tilde);
        glm::row(qmat, 1, ray2.q_hat);
        qmat = RotMat * qmat;
        vec2 q_tilde_in = glm::row(qmat, 0);
        vec2 q_hat_in   = glm::row(qmat, 1);
        
        // here's the actual curvature change
        
        vec2 p_tilde_out = p_tilde_in - q_tilde_in * r1 - q_hat_in * r2;
        vec2 p_hat_out   = p_hat_in   - q_tilde_in * r2;
        
        // rotate p back to e1, e2 system, q does not change
        // Note RotMat^(-1) = RotMat^T
        
        glm::row(pmat, 0, p_tilde_out);
        glm::row(pmat, 1, p_hat_out);
        pmat = glm::transpose(RotMat) * pmat;
        ray2.p_tilde = glm::row(pmat, 0);
        ray2.p_hat   = glm::row(pmat, 1);
        
    }else{
        // LP: 2D-3D only:
        // mbp: this needs modifying like the full 3D version to handle jumps in the x-y direction
        vec2 ray2n = vec2(-ray2.t.y, ray2.t.x); // ray normal
        
        real cnjump = glm::dot(gradcjump, ray2n);
        real csjump = glm::dot(gradcjump, ray2.t);
        
        real rm, rn;
        if(iSeg.z != iSeg0.z){         // crossing in depth
            rm =  ray2.t.x / ray2.t.y; // this is tan( alpha ) where alpha is the angle of incidence
        }else{                         // crossing in range
            rm = -ray2.t.y / ray2.t.x; // this is tan( alpha ) where alpha is the angle of incidence
        }                              // LP: The case where it crosses in depth and range simultaneously is not handled.
        
        rn = rm * (FL(2.0) * cnjump - rm * csjump) / ray2.c;
        ray2.p = ray2.p - ray2.q * rn;
    }
}

/**
 * Does a single step along the ray
 */
template<bool O3D, bool R3D> HOST_DEVICE inline void Step(
    rayPt<R3D> ray0, rayPt<R3D> &ray2, const Origin<O3D, R3D> &org,
    const BdryState<O3D> &bds, const BeamStructure *Beam, const SSPStructure *ssp,
    SSPSegState &iSeg, int32_t &iSmallStepCtr, bool &topRefl, bool &botRefl)
{
    rayPt<R3D> ray1;
    SSPOutputs<O3D> o0, o1, o2;
    StepPartials<R3D> part0, part1;
    rayPtExtras<R3D> pq0, pq1;
    VEC23<R3D> urayt0, urayt1;
    real csq0, csq1, h, w0, w1, hw0, hw1;
    
    #ifdef STEP_DEBUGGING
    printf("\nray0 x t p q tau amp (%20.17f,%20.17f) (%20.17f,%20.17f) (%20.17f,%20.17f) (%20.17f,%20.17f) (%20.17f,%20.17f) %20.17f\n", 
        ray0.x.x, ray0.x.y, ray0.t.x, ray0.t.y, ray0.p.x, ray0.p.y, ray0.q.x, ray0.q.y, ray0.tau.real(), ray0.tau.imag(), ray0.Amp);
    // printf("iSeg.z iSeg.r %d %d\n", iSeg.z, iSeg.r);
    #endif
    
    /*
    if(ray0.x.x > 16800.0){
        printf("Enough\n");
        bail();
    }
    */
    
    // The numerical integrator used here is a version of the polygon (a.k.a. midpoint, leapfrog, or Box method), and similar
    // to the Heun (second order Runge-Kutta method).
    // However, it's modified to allow for a dynamic step change, while preserving the second-order accuracy).

    // *** Phase 1 (an Euler step)

    EvaluateSSP<O3D, R3D>(ray0.x, ray0.t, o0, org, ssp, iSeg);
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
    ReduceStep<O3D>(x_o, t_o, iSeg0, bds, Beam, ssp, h, iSmallStepCtr);
    // printf("out h, urayt0 %20.17f (%20.17f, %20.17f)\n", h, urayt0.x, urayt0.y);
    real halfh = FL(0.5) * h; // first step of the modified polygon method is a half step
    
    ray1.x = ray0.x + halfh * urayt0;
    ray1.t = ray0.t - halfh * o0.gradc / csq0;
    UpdateRayPQ<R3D>(ray1, ray0, halfh, pq0);
    
    // printf("ray1 x t p q (%20.17f,%20.17f) (%20.17f,%20.17f) (%20.17f,%20.17f) (%20.17f,%20.17f)\n", 
    //     ray1.x.x, ray1.x.y, ray1.t.x, ray1.t.y, ray1.p.x, ray1.p.y, ray1.q.x, ray1.q.y);
    
    // *** Phase 2
    
    EvaluateSSP<O3D, R3D>(ray1.x, ray1.t, o1, ssp, iSeg);
    Get_c_partials<R3D>(ray1, o1, part1);
    pq1 = ComputeDeltaPQ<R3D>(ray1, o1, part1);
    
    // The Munk test case with a horizontally launched ray caused problems.
    // The ray vertexes on an interface and can ping-pong around that interface.
    // Have to be careful in that case about big changes to the stepsize (that invalidate the leap-frog scheme) in phase II.
    // A modified Heun or Box method could also work.

    csq1      = SQ(o1.ccpx.real());
    urayt1 = o1.ccpx.real() * ray1.t; // unit tangent
    
    // printf("urayt1 (%g,%g)\n", urayt1.x, urayt1.y);
    
    // reduce h to land on boundary
    t_o = RayToOceanT(urayt1, org);
    ReduceStep<O3D>(x_o, t_o, iSeg0, bds, Beam, ssp, h, iSmallStepCtr);
    
    // use blend of f' based on proportion of a full step used.
    w1  = h / (RL(2.0) * halfh);
    w0  = RL(1.0) - w1;
    // printf("w1 %20.17f w0 %20.17f\n", w1, w0);
    VEC23<R3D> urayt2 =  w0 * urayt0 + w1 * urayt1;
    // Take the blended ray tangent (urayt2) and find the minimum step size (h)
    // to put this on a boundary, and ensure that the resulting position
    // (ray2.x) gets put precisely on the boundary.
    VEC23<O3D> x2_o;
    StepToBdry<O3D>(x_o, x2_o, urayt2, h, topRefl, botRefl, iSeg0, bds, Beam, ssp);
    ray2.x = OceanToRayX(x2_o, org, urayt2);
    
    // Update other variables with this new h
    hw0 = h * w0;
    hw1 = h * w1;
    ray2.t   = ray0.t   - hw0 * o0.gradc / csq0         - hw1 * o1.gradc / csq1;
    ray2.tau = ray0.tau + hw0 / o0.ccpx                 + hw1 / o1.ccpx;
    UpdateRayPQ<R3D>(ray2, ray0, hw0, pq0);
    UpdateRayPQ<R3D>(ray2, ray2, hw1, pq1);
    
    // printf("ray2 x t p q tau (%g,%g) (%g,%g) (%g,%g) (%g,%g) (%g,%g)\n", 
    //     ray2.x.x, ray2.x.y, ray2.t.x, ray2.t.y, ray2.p.x, ray2.p.y, 
    //     ray2.q.x, ray2.q.y, ray2.tau.real(), ray2.tau.imag());
    
    ray2.Amp       = ray0.Amp;
    ray2.Phase     = ray0.Phase;
    ray2.NumTopBnc = ray0.NumTopBnc;
    ray2.NumBotBnc = ray0.NumBotBnc;
    
    // If we crossed an interface, apply jump condition

    EvaluateSSP<O3D, R3D>(ray2.x, ray2.t, o2, ssp, iSeg);
    ray2.c = o2.ccpx.real();
    
    if(iSeg.z != iSeg0.z || (!R3D && iSeg.r != iSeg0.r) || 
            (R3D && (iSeg.x != iSeg0.x || iSeg.y != iSeg0.y))){
        VEC23<R3D> gradcjump = o2.gradc - o0.gradc;
        CurvatureCorrection<R3D>(ray2, gradcjump, iSeg, iSeg0);
    }
}

}
