#pragma once
#include "common.hpp"

#define INFINITESIMAL_STEP_SIZE (RL(1e-6))

/**
 * calculate a reduced step size, h, that lands on any points where the environment changes
 * 
 * x0, urayt: ray coordinate and tangent
 * iSegz0, iSegr0: SSP layer the ray is in
 * Topx, Topn, Botx, Botn: Top, bottom coordinate and normal
 * h: reduced step size
 */
HOST_DEVICE inline void ReduceStep2D(const vec2 &x0, const vec2 &urayt,
    int32_t iSegz0, int32_t iSegr0, const vec2 &Topx, const vec2 &Topn,
    const vec2 &Botx, const vec2 &Botn, const vec2 &rTopSeg, const vec2 &rBotSeg,
    const BeamStructure *Beam, const SSPStructure *ssp, 
    real &h, int32_t &iSmallStepCtr)
{
    vec2 x, d, d0, rSeg;
    real h1, h2, h3, h4;
    
    // Detect interface or boundary crossing and reduce step, if necessary, to land on that crossing.
    // Keep in mind possibility that user put source right on an interface
    // and that multiple events can occur (crossing interface, top, and bottom in a single step).

    x = x0 + h * urayt; // make a trial step

    // interface crossing in depth
    h1 = REAL_MAX;
    if(STD::abs(urayt.y) > REAL_EPSILON){
        if(       ssp->z[iSegz0]     >  x.y){
            h1 = (ssp->z[iSegz0]     - x0.y) / urayt.y;
            // printf("Lower bound SSP Z %g > z %g; h1 = %g\n", ssp->z[iSegz0], x.y, h1);
        }else if( ssp->z[iSegz0 + 1] <  x.y){
            h1 = (ssp->z[iSegz0 + 1] - x0.y) / urayt.y;
            // printf("Upper bound SSP Z %g < z %g; h1 = %g\n", ssp->z[iSegz0+1], x.y, h1);
        }
    }
    
    // top crossing
    h2 = REAL_MAX;
    d = x - Topx; // vector from top to ray
    if(glm::dot(Topn, d) > REAL_EPSILON){
        d0 = x0 - Topx; // vector from top node to ray origin
        h2 = -glm::dot(d0, Topn) / glm::dot(urayt, Topn);
    }
    
    // bottom crossing
    h3 = REAL_MAX;
    d = x - Botx; // vector from bottom to ray
    if(glm::dot(Botn, d) > REAL_EPSILON){
        d0 = x0 - Botx; // vector from bottom node to ray origin
        h3 = -glm::dot(d0, Botn) / glm::dot(urayt, Botn);
    }
    
    // top or bottom segment crossing in range
    rSeg.x = math::max(rTopSeg.x, rBotSeg.x);
    rSeg.y = math::min(rTopSeg.y, rBotSeg.y);
    
    if(ssp->Type == 'Q'){
        rSeg.x = math::max(rSeg.x, ssp->Seg.r[iSegr0  ]);
        rSeg.y = math::min(rSeg.y, ssp->Seg.r[iSegr0+1]);
    }
    
    h4 = REAL_MAX;
    if(STD::abs(urayt.x) > REAL_EPSILON){
        if(x.x < rSeg.x){
            h4 = -(x0.x - rSeg.x) / urayt.x;
            // printf("Lower bound SSP R %g > r %g; h4 = %g\n", rSeg.x, x.x, h4);
        }else if(x.x > rSeg.y){
            h4 = -(x0.x - rSeg.y) / urayt.x;
            // printf("Upper bound SSP R %g < r %g; h4 = %g; rTopSeg up %g; rBotSeg up %g; ssp r up %g\n",
            //     rSeg.y, x.x, h4, rTopSeg.y, rBotSeg.y,
            //     ssp->Type == 'Q' ? ssp->Seg.r[iSegr0+1] : INFINITY);
        }
    }
    
    //printf("ReduceStep2D h h1 h2 h3 h4 %g %g %g %g %g\n", h, h1, h2, h3, h4);
    // take limit set by shortest distance to a crossing
    h = math::min(h, math::min(h1, math::min(h2, math::min(h3, h4))));
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
    int32_t iSegz0, int32_t iSegr0, const vec2 &Topx, const vec2 &Topn,
    const vec2 &Botx, const vec2 &Botn, const vec2 &rTopSeg, const vec2 &rBotSeg,
    const BeamStructure *Beam, const SSPStructure *ssp)
{
    vec2 d, d0, rSeg;
    
    // Original step due to maximum step size
    h = Beam->deltas;
    x2 = x0 + h * urayt;
    
    // interface crossing in depth
    if(STD::abs(urayt.y) > REAL_EPSILON){
        if(      ssp->z[iSegz0]     > x2.y){
            h = (ssp->z[iSegz0]     - x0.y) / urayt.y;
            x2.x = x0.x + h * urayt.x;
            x2.y = ssp->z[iSegz0];
            // printf("StepToBdry2D lower depth h %g to (%g,%g)\n", h, x2.x, x2.y);
        }else if(ssp->z[iSegz0 + 1] < x2.y){
            h = (ssp->z[iSegz0 + 1] - x0.y) / urayt.y;
            x2.x = x0.x + h * urayt.x;
            x2.y = ssp->z[iSegz0 + 1];
            // printf("StepToBdry2D upper depth h %g to (%g,%g)\n", h, x2.x, x2.y);
        }
    }
    
    // top crossing
    d = x2 - Topx; // vector from top to ray
    if(glm::dot(Topn, d) > REAL_EPSILON){
        d0 = x0 - Topx; // vector from top node to ray origin
        h = -glm::dot(d0, Topn) / glm::dot(urayt, Topn);
        x2 = x0 + h * urayt;
        // Snap to exact top depth value if it's flat
        if(STD::abs(Topn.x) < REAL_EPSILON){
            x2.y = Topx.y;
        }
        // printf("StepToBdry2D top crossing h %g to (%g,%g)\n", h, x2.x, x2.y);
        topRefl = true;
    }else{
        topRefl = false;
    }
    
    // bottom crossing
    d = x2 - Botx; // vector from bottom to ray
    if(glm::dot(Botn, d) > REAL_EPSILON){
        d0 = x0 - Botx; // vector from bottom node to ray origin
        h = -glm::dot(d0, Botn) / glm::dot(urayt, Botn);
        x2 = x0 + h * urayt;
        // Snap to exact bottom depth value if it's flat
        if(STD::abs(Botn.x) < REAL_EPSILON){
            x2.y = Botx.y;
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
    rSeg.x = math::max(rTopSeg.x, rBotSeg.x); // lower range bound (not an x value)
    rSeg.y = math::min(rTopSeg.y, rBotSeg.y); // upper range bound (not a y value)
    
    if(ssp->Type == 'Q'){
        rSeg.x = math::max(rSeg.x, ssp->Seg.r[iSegr0  ]);
        rSeg.y = math::min(rSeg.y, ssp->Seg.r[iSegr0+1]);
    }
    
    if(STD::abs(urayt.x) > REAL_EPSILON){
        if(x2.x < rSeg.x){
            h = -(x0.x - rSeg.x) / urayt.x;
            x2.x = rSeg.x;
            x2.y = x0.y + h * urayt.y;
            // printf("StepToBdry2D lower range h %g to (%g,%g)\n", h, x2.x, x2.y);
            topRefl = false;
            botRefl = false;
        }else if(x2.x > rSeg.y){
            h = -(x0.x - rSeg.y) / urayt.x;
            x2.x = rSeg.y;
            x2.y = x0.y + h * urayt.y;
            // printf("StepToBdry2D upper range h %g to (%g,%g)\n", h, x2.x, x2.y);
            topRefl = false;
            botRefl = false;
        }
    }
    
    if(h < INFINITESIMAL_STEP_SIZE * Beam->deltas){ // is it taking an infinitesimal step?
        h = INFINITESIMAL_STEP_SIZE * Beam->deltas; // make sure we make some motion
        x2 = x0 + h * urayt;
        // printf("StepToBdry2D small step forced h %g to (%g,%g)\n", h, x2.x, x2.y);
        // Recheck reflection conditions
        d = x2 - Topx; // vector from top to ray
        if(glm::dot(Topn, d) > REAL_EPSILON){
            topRefl = true;
        }else{
            topRefl = false;
        }
        d = x2 - Botx; // vector from bottom to ray
        if(glm::dot(Botn, d) > REAL_EPSILON){
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
    const vec2 &Topx, const vec2 &Topn, const vec2 &Botx, const vec2 &Botn,
    const vec2 &rTopSeg, const vec2 &rBotSeg, const real &freq,
    const BeamStructure *Beam, const SSPStructure *ssp,
    int32_t &iSegz, int32_t &iSegr, int32_t &iSmallStepCtr, bool &topRefl, bool &botRefl)
{
    ray2DPt ray1;
    int32_t iSegz0, iSegr0;
    vec2 gradc0, gradc1, gradc2, urayt0, urayt1, ray2n, gradcjump;
    cpx ccpx0, ccpx1, ccpx2;
    real crr0, crz0, czz0, csq0, cnn0_csq0;
    real crr1, crz1, czz1, csq1, cnn1_csq1;
    real crr2, crz2, czz2;
    real h, halfh, rm, rn, cnjump, csjump, w0, w1, rho;
    
    // printf("\nray0 x t p q tau amp (%20.17f,%20.17f) (%20.17f,%20.17f) (%20.17f,%20.17f) (%20.17f,%20.17f) (%20.17f,%20.17f) %20.17f\n", 
    //     ray0.x.x, ray0.x.y, ray0.t.x, ray0.t.y, ray0.p.x, ray0.p.y, ray0.q.x, ray0.q.y, ray0.tau.real(), ray0.tau.imag(), ray0.Amp);
    // printf("iSegz iSegr %d %d\n", iSegz, iSegr);
    
    // if(ray0.x.x > 420.0){
    //     printf("Enough\n");
    //     bail();
    // }
    
    // The numerical integrator used here is a version of the polygon (a.k.a. midpoint, leapfrog, or Box method), and similar
    // to the Heun (second order Runge-Kutta method).
    // However, it's modified to allow for a dynamic step change, while preserving the second-order accuracy).

    // *** Phase 1 (an Euler step)

    EvaluateSSP(ray0.x, ray0.t, ccpx0, gradc0, crr0, crz0, czz0, rho, freq, ssp, iSegz, iSegr);
    // printf("iSegz iSegr %d %d\n", iSegz, iSegr);
    
    // printf("ccpx0: (%g,%g) gradc0: (%g,%g) crr0 %g crz0 %g czz0 %g, rho %g\n", 
    //    ccpx0.real(), ccpx0.imag(), gradc0.x, gradc0.y, crr0, crz0, czz0, rho);
    
    /*
    if(STD::abs(ccpx0) > DEBUG_LARGEVAL){
        printf("ccpx0 invalid: (%g,%g)\n", ccpx0.real(), ccpx0.imag());
        bail();
    }
    */
    
    csq0      = SQ(ccpx0.real());
    cnn0_csq0 = crr0 * SQ(ray0.t.y) - FL(2.0) * crz0 * ray0.t.x * ray0.t.y + czz0 * SQ(ray0.t.x);
    iSegz0    = iSegz; // make note of current layer
    iSegr0    = iSegr;
    
    h = Beam->deltas;       // initially set the step h, to the basic one, deltas
    urayt0 = ccpx0.real() * ray0.t; // unit tangent
    
    // printf("urayt0 (%g,%g)\n", urayt0.x, urayt0.y);
    
    // reduce h to land on boundary
    ReduceStep2D(ray0.x, urayt0, iSegz0, iSegr0, Topx, Topn, Botx, Botn, rTopSeg, rBotSeg,
        Beam, ssp, h, iSmallStepCtr);
    // printf("out h, urayt0 %20.17f (%20.17f, %20.17f)\n", h, urayt0.x, urayt0.y);
    halfh = FL(0.5) * h; // first step of the modified polygon method is a half step
    
    ray1.x = ray0.x + halfh * urayt0;
    ray1.t = ray0.t - halfh * gradc0 / csq0;
    ray1.p = ray0.p - halfh * cnn0_csq0    * ray0.q;
    ray1.q = ray0.q + halfh * ccpx0.real() * ray0.p;
    
    // printf("ray1 x t p q (%20.17f,%20.17f) (%20.17f,%20.17f) (%20.17f,%20.17f) (%20.17f,%20.17f)\n", 
    //   ray1.x.x, ray1.x.y, ray1.t.x, ray1.t.y, ray1.p.x, ray1.p.y, ray1.q.x, ray1.q.y);
    
    /*
    if(STD::abs(ray1.x.x) > DEBUG_LARGEVAL || STD::abs(ray1.x.y) > DEBUG_LARGEVAL){
        printf("ray1.x invalid\n");
        bail();
    }
    if(STD::abs(ray1.t.x) > DEBUG_LARGEVAL || STD::abs(ray1.t.y) > DEBUG_LARGEVAL){
        printf("ray1.t invalid\n");
        bail();
    }
    if(STD::abs(ray1.p.x) > DEBUG_LARGEVAL || STD::abs(ray1.p.y) > DEBUG_LARGEVAL){
        printf("ray1.p invalid\n");
        bail();
    }
    if(STD::abs(ray1.q.x) > DEBUG_LARGEVAL || STD::abs(ray1.q.y) > DEBUG_LARGEVAL){
        printf("ray1.q invalid\n");
        bail();
    }
    */
    
    // *** Phase 2
    
    EvaluateSSP(ray1.x, ray1.t, ccpx1, gradc1, crr1, crz1, czz1, rho, freq, ssp, iSegz, iSegr);
    
    // printf("ccpx1: (%g,%g)\n", ccpx1.real(), ccpx1.imag());
    
    /*
    if(STD::abs(ccpx1) > DEBUG_LARGEVAL){
        printf("ccpx1 invalid: ray1.x (%g,%g) iSegz %d iSegr %d => ccpx1 (%g,%g)\n", 
            ray1.x.x, ray1.x.y, iSegz, iSegr, ccpx1.real(), ccpx1.imag());
        bail();
    }
    */
    csq1      = SQ(ccpx1.real());
    cnn1_csq1 = crr1 * SQ(ray1.t.y) - FL(2.0) * crz1 * ray1.t.x * ray1.t.y + czz1 * SQ(ray1.t.x);
    
    // The Munk test case with a horizontally launched ray caused problems.
    // The ray vertexes on an interface and can ping-pong around that interface.
    // Have to be careful in that case about big changes to the stepsize (that invalidate the leap-frog scheme) in phase II.
    // A modified Heun or Box method could also work.

    urayt1 = ccpx1.real() * ray1.t; // unit tangent
    
    // printf("urayt1 (%g,%g)\n", urayt1.x, urayt1.y);
    
    // reduce h to land on boundary
    ReduceStep2D(ray0.x, urayt1, iSegz0, iSegr0, Topx, Topn, Botx, Botn, rTopSeg, rBotSeg,
        Beam, ssp, h, iSmallStepCtr);
    
    // use blend of f' based on proportion of a full step used.
    w1  = h / (RL(2.0) * halfh);
    w0  = RL(1.0) - w1;
    // printf("w1 %20.17f w0 %20.17f\n", w1, w0);
    vec2 urayt2   =  w0 * urayt0                + w1 * urayt1;
    vec2 unitdt   = -w0 * gradc0 / csq0         - w1 * gradc1 / csq1;
    vec2 unitdp   = -w0 * cnn0_csq0    * ray0.q - w1 * cnn1_csq1    * ray1.q;
    vec2 unitdq   =  w0 * ccpx0.real() * ray0.p + w1 * ccpx1.real() * ray1.p;
    cpx  unitdtau =  w0 / ccpx0                 + w1 / ccpx1;
    
    // Take the blended ray tangent (urayt2) and find the minimum step size (h)
    // to put this on a boundary, and ensure that the resulting position
    // (ray2.x) gets put precisely on the boundary.
    StepToBdry2D(ray0.x, ray2.x, urayt2, h, topRefl, botRefl,
        iSegz0, iSegr0, Topx, Topn, Botx, Botn, rTopSeg, rBotSeg, Beam, ssp);
    ray2.t   = ray0.t   + h * unitdt;
    ray2.p   = ray0.p   + h * unitdp;
    ray2.q   = ray0.q   + h * unitdq;
    ray2.tau = ray0.tau + h * unitdtau;
    
    /*
    printf("ray2 x t p q tau (%g,%g) (%g,%g) (%g,%g) (%g,%g) (%g,%g)\n", 
        ray2.x.x, ray2.x.y, ray2.t.x, ray2.t.y, ray2.p.x, ray2.p.y, 
        ray2.q.x, ray2.q.y, ray2.tau.real(), ray2.tau.imag());
    
    if(STD::abs(ray2.x.x) > DEBUG_LARGEVAL || STD::abs(ray2.x.y) > DEBUG_LARGEVAL
    || STD::abs(ray2.t.x) > DEBUG_LARGEVAL || STD::abs(ray2.t.y) > DEBUG_LARGEVAL
    || STD::abs(ray2.p.x) > DEBUG_LARGEVAL || STD::abs(ray2.p.y) > DEBUG_LARGEVAL
    || STD::abs(ray2.q.x) > DEBUG_LARGEVAL || STD::abs(ray2.q.y) > DEBUG_LARGEVAL){
        printf("ray2 invalid\n");
        bail();
    }
    */
    
    ray2.Amp       = ray0.Amp;
    ray2.Phase     = ray0.Phase;
    ray2.NumTopBnc = ray0.NumTopBnc;
    ray2.NumBotBnc = ray0.NumBotBnc;
    
    // If we crossed an interface, apply jump condition

    EvaluateSSP(ray2.x, ray2.t, ccpx2, gradc2, crr2, crz2, czz2, rho, freq, ssp, iSegz, iSegr);
    //printf("ccpx2: (%g,%g)\n", ccpx2.real(), ccpx2.imag());
    /*
    if(STD::abs(ccpx2) > DEBUG_LARGEVAL){
        printf("ccpx2 invalid: (%g,%g)\n", ccpx1.real(), ccpx1.imag());
        bail();
    }
    */
    ray2.c = ccpx2.real();
    
    if(iSegz != iSegz0 || iSegr != iSegr0){
        gradcjump = gradc2 - gradc0;
        ray2n = vec2(-ray2.t.y, ray2.t.x); // ray normal
        
        cnjump = glm::dot(gradcjump, ray2n);
        csjump = glm::dot(gradcjump, ray2.t);
        
        if(iSegz != iSegz0){             // crossing in depth
            rm =  ray2.t.x / ray2.t.y; // this is tan( alpha ) where alpha is the angle of incidence
        }else{                           // crossing in range
            rm = -ray2.t.y / ray2.t.x; // this is tan( alpha ) where alpha is the angle of incidence
        }                                // LP: The case where it crosses in depth and range simultaneously is not handled.
        
        rn = rm * (FL(2.0) * cnjump - rm * csjump) / ccpx2.real();
        ray2.p = ray2.p - ray2.q * rn;
    }
}
