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

#ifdef BHC_USE_FLOATS
#define INFINITESIMAL_STEP_SIZE (RL(1e-3))
#else
#define INFINITESIMAL_STEP_SIZE (RL(1e-6))
#endif

/**
 * top or bottom segment crossing in range / x / y
 */
HOST_DEVICE inline real TopBotSegCrossing(
    const BdryLimits &TopSeg, const BdryLimits &BotSeg,
    const real *segdim, int32_t iSeg,
    real xdim, real x0dim, real uraytdim, char qh, const SSPStructure *ssp
    )
{
    BdryLimits rSeg;
    rSeg.min = bhc::max(TopSeg.min, BotSeg.min);
    rSeg.max = bhc::min(TopSeg.max, BotSeg.max);
    
    if(ssp->Type == qh){
        rSeg.min = bhc::max(rSeg.min, segdim[iSeg  ]);
        rSeg.max = bhc::min(rSeg.max, segdim[iSeg+1]);
    }
    
    real h = REAL_MAX;
    if(STD::abs(uraytdim) > REAL_EPSILON){
        if(xdim < rSeg.min){
            h = -(x0dim - rSeg.min) / uraytdim;
            // printf("Closer bound SSP R %g > r %g; h4 = %g\n", rSeg.min, xdim, h);
        }else if(x.x > rSeg.max){
            h = -(x0dim - rSeg.max) / uraytdim;
            // printf("Farther bound SSP R %g < r %g; h4 = %g; top.max %g; bot.max %g; ssp r up %g\n",
            //     rSeg.max, xdim, h, TopSeg.max, BotSeg.max,
            //     ssp->Type == qh ? segdim[iSeg+1] : INFINITY);
        }
    }
    
    return h;
}

/**
 * calculate a reduced step size, h, that lands on any points where the environment changes
 * 
 * x0, urayt: ray coordinate and tangent
 * iSeg0: SSP layer the ray is in
 * Topx, Topn, Botx, Botn: Top, bottom coordinate and normal
 * h: reduced step size
 */
template<bool THREED> HOST_DEVICE inline void ReduceStep2D(
    const typename TmplVec23<THREED>::type &x0, 
    const typename TmplVec23<THREED>::type &urayt,
    const SSPSegState &iSeg0, const BdryState<THREED> &bds,
    const BeamStructure *Beam, const SSPStructure *ssp, 
    real &h, int32_t &iSmallStepCtr)
{
    typename TmplVec23<THREED>::type x, d, d0;
    real h1, h2, h3, h4, h5, h6, h7;
    
    // Detect interface or boundary crossing and reduce step, if necessary, to land on that crossing.
    // Keep in mind possibility that user put source right on an interface
    // and that multiple events can occur (crossing interface, top, and bottom in a single step).

    x = x0 + h * urayt; // make a trial step

    // interface crossing in depth
    // LP: 3D only:
    // Step reduction is not done for the top or bottom layer
    // Instead the SSP is extrapolated
    // This prevents problems when the boundaries are outside the domain of the SSP
    h1 = REAL_MAX;
    if(STD::abs(DEP(urayt)) > REAL_EPSILON){
        if(       ssp->z[iSeg0.z]     > DEP(x) && (!THREED || iSeg0.z > 0)){
            h1 = (ssp->z[iSeg0.z]     - DEP(x0)) / DEP(urayt);
            // printf("Shallower bound SSP Z %g > z %g; h1 = %g\n", ssp->z[iSeg0.z], x.y, h1);
        }else if( ssp->z[iSeg0.z + 1] < DEP(x) && (!THREED || iSeg0.z + 1 < ssp->Nz-1)){
            h1 = (ssp->z[iSeg0.z + 1] - DEP(x0)) / DEP(urayt);
            // printf("Deeper bound SSP Z %g < z %g; h1 = %g\n", ssp->z[iSeg0.z+1], x.y, h1);
        }
    }
    
    // top crossing
    h2 = REAL_MAX;
    d = x - bds.top.x; // vector from top to ray
    if(glm::dot(bds.top.n, d) >= RL(0.0)){
        d0 = x0 - bds.top.x; // vector from top node to ray origin
        h2 = -glm::dot(d0, bds.top.n) / glm::dot(urayt, bds.top.n);
    }
    
    // bottom crossing
    h3 = REAL_MAX;
    d = x - bds.bot.x; // vector from bottom to ray
    if(glm::dot(bds.bot.n, d) >= RL(0.0)){
        d0 = x0 - bds.bot.x; // vector from bottom node to ray origin
        h3 = -glm::dot(d0, bds.bot.n) / glm::dot(urayt, bds.bot.n);
    }
    
    if constexpr(THREED){
        h4 = TopBotSegCrossing(bds.top.lSeg.x, bds.bot.lSeg.x, ssp->Seg.x, iSeg0.x,
            x.x, x0.x, urayt.x, 'H', ssp);
        h5 = TopBotSegCrossing(bds.top.lSeg.y, bds.bot.lSeg.y, ssp->Seg.y, iSeg0.y,
            x.y, x0.y, urayt.y, 'H', ssp);
    }else{
        h4 = TopBotSegCrossing(bds.top.lSeg, bds.bot.lSeg, ssp->Seg.r, iSeg0.r,
            x.x, x0.x, urayt.x, 'Q', ssp);
        h5 = REAL_MAX;
    }
    
    //printf("ReduceStep2D h h1 h2 h3 h4 h5 h6 h7 %g %g %g %g %g %g %g %g\n",
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
        printf("Step %g due to R SSP crossing\n", h);
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

HOST_DEVICE inline void StepToBdry2D(const vec2 &x0, vec2 &x2, const vec2 &urayt,
    real &h, bool &topRefl, bool &botRefl,
    const SSPSegState &iSeg0, const BdryState<false> &bds,
    const BeamStructure *Beam, const SSPStructure *ssp)
{
    vec2 d, d0;
    BdryLimits rSeg;
    
    // Original step due to maximum step size
    h = Beam->deltas;
    x2 = x0 + h * urayt;
    
    // interface crossing in depth
    if(STD::abs(urayt.y) > REAL_EPSILON){
        if(      ssp->z[iSeg0.z]     > x2.y){
            h = (ssp->z[iSeg0.z]     - x0.y) / urayt.y;
            x2.x = x0.x + h * urayt.x;
            x2.y = ssp->z[iSeg0.z];
            // printf("StepToBdry2D upper depth h %g to (%g,%g)\n", h, x2.x, x2.y);
        }else if(ssp->z[iSeg0.z + 1] < x2.y){
            h = (ssp->z[iSeg0.z + 1] - x0.y) / urayt.y;
            x2.x = x0.x + h * urayt.x;
            x2.y = ssp->z[iSeg0.z + 1];
            // printf("StepToBdry2D lower depth h %g to (%g,%g)\n", h, x2.x, x2.y);
        }
    }
    
    // top crossing
    d = x2 - bds.top.x; // vector from top to ray
    // Originally, this value had to be > a small positive number, meaning the
    // new step really had to be outside the boundary, not just to the boundary.
    // Also, this is not missing a normalization factor, Topn is normalized so
    // this is actually the distance above the top in meters.
    if(glm::dot(bds.top.n, d) > -INFINITESIMAL_STEP_SIZE){
        d0 = x0 - bds.top.x; // vector from top node to ray origin
        h = -glm::dot(d0, bds.top.n) / glm::dot(urayt, bds.top.n);
        x2 = x0 + h * urayt;
        // Snap to exact top depth value if it's flat
        if(STD::abs(bds.top.n.x) < REAL_EPSILON){
            x2.y = bds.top.x.y;
        }
        // printf("StepToBdry2D top crossing h %g to (%g,%g)\n", h, x2.x, x2.y);
        topRefl = true;
    }else{
        topRefl = false;
    }
    
    // bottom crossing
    d = x2 - bds.bot.x; // vector from bottom to ray
    // See comment above for top case.
    if(glm::dot(bds.bot.n, d) > -INFINITESIMAL_STEP_SIZE){
        d0 = x0 - bds.bot.x; // vector from bottom node to ray origin
        h = -glm::dot(d0, bds.bot.n) / glm::dot(urayt, bds.bot.n);
        x2 = x0 + h * urayt;
        // Snap to exact bottom depth value if it's flat
        if(STD::abs(bds.bot.n.x) < REAL_EPSILON){
            x2.y = bds.bot.x.y;
        }
        // printf("StepToBdry2D bottom crossing h %g to (%g,%g)\n", h, x2.x, x2.y);
        botRefl = true;
        // Should not ever be able to cross both, but in case it does, make sure
        // only the crossing we exactly landed on is active
        topRefl = false;
    }else{
        botRefl = false;
    }
    
    // top or bottom segment crossing in range
    rSeg.min = bhc::max(bds.top.lSeg.min, bds.bot.lSeg.min);
    rSeg.max = bhc::min(bds.top.lSeg.max, bds.bot.lSeg.max);
    
    if(ssp->Type == 'Q'){
        rSeg.min = bhc::max(rSeg.min, ssp->Seg.r[iSeg0.r  ]);
        rSeg.max = bhc::min(rSeg.max, ssp->Seg.r[iSeg0.r+1]);
    }
    
    if(STD::abs(urayt.x) > REAL_EPSILON){
        if(x2.x < rSeg.min){
            h = -(x0.x - rSeg.min) / urayt.x;
            x2.x = rSeg.min;
            x2.y = x0.y + h * urayt.y;
            // printf("StepToBdry2D lower range h %g to (%g,%g)\n", h, x2.x, x2.y);
            topRefl = false;
            botRefl = false;
        }else if(x2.x > rSeg.max){
            h = -(x0.x - rSeg.max) / urayt.x;
            x2.x = rSeg.max;
            x2.y = x0.y + h * urayt.y;
            // printf("StepToBdry2D upper range h %25.21f to (%25.21f,%25.21f)\n", h, x2.x, x2.y);
            topRefl = false;
            botRefl = false;
        }
    }
    
    if(h < INFINITESIMAL_STEP_SIZE * Beam->deltas){ // is it taking an infinitesimal step?
        h = INFINITESIMAL_STEP_SIZE * Beam->deltas; // make sure we make some motion
        x2 = x0 + h * urayt;
        // printf("StepToBdry2D small step forced h %g to (%g,%g)\n", h, x2.x, x2.y);
        // Recheck reflection conditions
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

/**
 * Does a single step along the ray
 */
HOST_DEVICE inline void Step2D(ray2DPt ray0, ray2DPt &ray2, 
    const BdryState<false> &bds, const BeamStructure *Beam, const SSPStructure *ssp,
    SSPSegState &iSeg, int32_t &iSmallStepCtr, bool &topRefl, bool &botRefl)
{
    ray2DPt ray1;
    SSPSegState iSeg0;
    vec2 urayt0, urayt1, ray2n, gradcjump;
    real csq0, cnn0_csq0, csq1, cnn1_csq1;
    real h, halfh, rm, rn, cnjump, csjump, w0, w1;
    SSPOutputs<false> o0, o1, o2;
    
    // printf("\nray0 x t p q tau amp (%20.17f,%20.17f) (%20.17f,%20.17f) (%20.17f,%20.17f) (%20.17f,%20.17f) (%20.17f,%20.17f) %20.17f\n", 
    //     ray0.x.x, ray0.x.y, ray0.t.x, ray0.t.y, ray0.p.x, ray0.p.y, ray0.q.x, ray0.q.y, ray0.tau.real(), ray0.tau.imag(), ray0.Amp);
    // printf("iSeg.z iSeg.r %d %d\n", iSeg.z, iSeg.r);
    
    // if(ray0.x.x > 40.0){
    //     printf("Enough\n");
    //     bail();
    // }
    
    // The numerical integrator used here is a version of the polygon (a.k.a. midpoint, leapfrog, or Box method), and similar
    // to the Heun (second order Runge-Kutta method).
    // However, it's modified to allow for a dynamic step change, while preserving the second-order accuracy).

    // *** Phase 1 (an Euler step)

    EvaluateSSP<false>(ray0.x, ray0.t, o0, ssp, iSeg);
    // printf("iSeg.z iSeg.r %d %d\n", iSeg.z, iSeg.r);
    
    csq0      = SQ(o0.ccpx.real());
    cnn0_csq0 = o0.crr * SQ(ray0.t.y) - FL(2.0) * o0.crz * ray0.t.x * ray0.t.y + o0.czz * SQ(ray0.t.x);
    iSeg0     = iSeg; // make note of current layer
    
    h = Beam->deltas;       // initially set the step h, to the basic one, deltas
    urayt0 = o0.ccpx.real() * ray0.t; // unit tangent
    
    // printf("urayt0 (%g,%g)\n", urayt0.x, urayt0.y);
    
    // reduce h to land on boundary
    ReduceStep2D(ray0.x, urayt0, iSeg0, bds, Beam, ssp, h, iSmallStepCtr);
    // printf("out h, urayt0 %20.17f (%20.17f, %20.17f)\n", h, urayt0.x, urayt0.y);
    halfh = FL(0.5) * h; // first step of the modified polygon method is a half step
    
    ray1.x = ray0.x + halfh * urayt0;
    ray1.t = ray0.t - halfh * o0.gradc / csq0;
    ray1.p = ray0.p - halfh * cnn0_csq0      * ray0.q;
    ray1.q = ray0.q + halfh * o0.ccpx.real() * ray0.p;
    
    // printf("ray1 x t p q (%20.17f,%20.17f) (%20.17f,%20.17f) (%20.17f,%20.17f) (%20.17f,%20.17f)\n", 
    //   ray1.x.x, ray1.x.y, ray1.t.x, ray1.t.y, ray1.p.x, ray1.p.y, ray1.q.x, ray1.q.y);
    
    // *** Phase 2
    
    EvaluateSSP<false>(ray1.x, ray1.t, o1, ssp, iSeg);
    
    csq1      = SQ(o1.ccpx.real());
    cnn1_csq1 = o1.crr * SQ(ray1.t.y) - FL(2.0) * o1.crz * ray1.t.x * ray1.t.y + o1.czz * SQ(ray1.t.x);
    
    // The Munk test case with a horizontally launched ray caused problems.
    // The ray vertexes on an interface and can ping-pong around that interface.
    // Have to be careful in that case about big changes to the stepsize (that invalidate the leap-frog scheme) in phase II.
    // A modified Heun or Box method could also work.

    urayt1 = o1.ccpx.real() * ray1.t; // unit tangent
    
    // printf("urayt1 (%g,%g)\n", urayt1.x, urayt1.y);
    
    // reduce h to land on boundary
    ReduceStep2D(ray0.x, urayt1, iSeg0, bds, Beam, ssp, h, iSmallStepCtr);
    
    // use blend of f' based on proportion of a full step used.
    w1  = h / (RL(2.0) * halfh);
    w0  = RL(1.0) - w1;
    // printf("w1 %20.17f w0 %20.17f\n", w1, w0);
    vec2 urayt2   =  w0 * urayt0                  + w1 * urayt1;
    vec2 unitdt   = -w0 * o0.gradc / csq0         - w1 * o1.gradc / csq1;
    vec2 unitdp   = -w0 * cnn0_csq0      * ray0.q - w1 * cnn1_csq1      * ray1.q;
    vec2 unitdq   =  w0 * o0.ccpx.real() * ray0.p + w1 * o1.ccpx.real() * ray1.p;
    cpx  unitdtau =  w0 / o0.ccpx                 + w1 / o1.ccpx;
    
    // Take the blended ray tangent (urayt2) and find the minimum step size (h)
    // to put this on a boundary, and ensure that the resulting position
    // (ray2.x) gets put precisely on the boundary.
    StepToBdry2D(ray0.x, ray2.x, urayt2, h, topRefl, botRefl, iSeg0, bds, Beam, ssp);
    ray2.t   = ray0.t   + h * unitdt;
    ray2.p   = ray0.p   + h * unitdp;
    ray2.q   = ray0.q   + h * unitdq;
    ray2.tau = ray0.tau + h * unitdtau;
    
    /*
    printf("ray2 x t p q tau (%g,%g) (%g,%g) (%g,%g) (%g,%g) (%g,%g)\n", 
        ray2.x.x, ray2.x.y, ray2.t.x, ray2.t.y, ray2.p.x, ray2.p.y, 
        ray2.q.x, ray2.q.y, ray2.tau.real(), ray2.tau.imag());
    */
    
    ray2.Amp       = ray0.Amp;
    ray2.Phase     = ray0.Phase;
    ray2.NumTopBnc = ray0.NumTopBnc;
    ray2.NumBotBnc = ray0.NumBotBnc;
    
    // If we crossed an interface, apply jump condition

    EvaluateSSP<false>(ray2.x, ray2.t, o2, ssp, iSeg);
    ray2.c = o2.ccpx.real();
    
    if(iSeg.z != iSeg0.z || iSeg.r != iSeg0.r){
        gradcjump = o2.gradc - o0.gradc;
        ray2n = vec2(-ray2.t.y, ray2.t.x); // ray normal
        
        cnjump = glm::dot(gradcjump, ray2n);
        csjump = glm::dot(gradcjump, ray2.t);
        
        if(iSeg.z != iSeg0.z){         // crossing in depth
            rm =  ray2.t.x / ray2.t.y; // this is tan( alpha ) where alpha is the angle of incidence
        }else{                         // crossing in range
            rm = -ray2.t.y / ray2.t.x; // this is tan( alpha ) where alpha is the angle of incidence
        }                              // LP: The case where it crosses in depth and range simultaneously is not handled.
        
        rn = rm * (FL(2.0) * cnjump - rm * csjump) / o2.ccpx.real();
        ray2.p = ray2.p - ray2.q * rn;
    }
}

}
