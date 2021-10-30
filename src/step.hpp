#pragma once
#include "common.hpp"
#include "ssp.hpp"
#include "beams.hpp"

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
    real &h, uint32_t &iSmallStepCtr)
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
        }else if( ssp->z[iSegz0 + 1] <  x.y){
            h1 = (ssp->z[iSegz0 + 1] - x0.y) / urayt.y;
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
    rSeg.x = STD::max(rTopSeg.x, rBotSeg.x);
    rSeg.y = STD::min(rTopSeg.y, rBotSeg.y);
    
    if(ssp->Type == 'Q'){
        rSeg.x = STD::max(rSeg.x, ssp->Seg.r[iSegr0  ]);
        rSeg.y = STD::min(rSeg.y, ssp->Seg.r[iSegr0+1]);
    }
    
    h4 = REAL_MAX;
    if(STD::abs(urayt.x) > REAL_EPSILON){
        if(x.x < rSeg.x){
            h4 = -(x0.x - rSeg.x) / urayt.x;
        }else if(x.x > rSeg.y){
            h4 = -(x0.x - rSeg.y) / urayt.x;
        }
    }
    
    // take limit set by shortest distance to a crossing
    h = STD::min(h, STD::min(h1, STD::min(h2, STD::min(h3, h4))));
    if(h < RC(1.0e-4) * Beam->deltas){ // is it taking an infinitesimal step?
        h = RC(1.0e-5) * Beam->deltas; // make sure we make some motion
        ++iSmallStepCtr; // keep a count of the number of sequential small steps
    }else{
        iSmallStepCtr = 0;
    }
}

/**
 * Does a single step along the ray
 */
HOST_DEVICE inline void Step2D(ray2DPt ray0, ray2DPt *ray2, 
    const vec2 &Topx, const vec2 &Topn, const vec2 &Botx, const vec2 &Botn,
    const vec2 &rTopSeg, const vec2 &rBotSeg, const real &freq,
    const BeamStructure *Beam, const SSPStructure *ssp,
    int32_t &iSegz, int32_t &iSegr, uint32_t &iSmallStepCtr)
{
    ray2DPt ray1;
    int32_t iSegz0, iSegr0;
    vec2 gradc0, gradc1, gradc2, urayt0, urayt1, ray2n, gradcjump;
    cpx ccpx0, ccpx1, ccpx2;
    real crr0, crz0, czz0, csq0, cnn0_csq0;
    real crr1, crz1, czz1, csq1, cnn1_csq1;
    real crr2, crz2, czz2;
    real h, halfh, hw0, hw1, rm, rn, cnjump, csjump, w0, w1, rho;
    
    // The numerical integrator used here is a version of the polygon (a.k.a. midpoint, leapfrog, or Box method), and similar
    // to the Heun (second order Runge-Kutta method).
    // However, it's modified to allow for a dynamic step change, while preserving the second-order accuracy).

    // *** Phase 1 (an Euler step)

    EvaluateSSP(ray0.x, ccpx0, gradc0, crr0, crz0, czz0, rho, freq, ssp, iSegz, iSegr);
    
    csq0      = SQ(ccpx0.real());
    cnn0_csq0 = crr0 * SQ(ray0.t.y) - RC(2.0) * crz0 * ray0.t.x * ray0.t.y + czz0 * SQ(ray0.t.x);
    iSegz0    = iSegz; // make note of current layer
    iSegr0    = iSegr;
    
    h = Beam->deltas;       // initially set the step h, to the basic one, deltas
    urayt0 = ccpx0.real() * ray0.t; // unit tangent
    
    // reduce h to land on boundary
    ReduceStep2D(ray0.x, urayt0, iSegz0, iSegr0, Topx, Topn, Botx, Botn, rTopSeg, rBotSeg,
        Beam, ssp, h, iSmallStepCtr);
    halfh = RC(0.5) * h; // first step of the modified polygon method is a half step
    
    ray1.x = ray0.x + halfh * urayt0;
    ray1.t = ray0.t - halfh * gradc0 / csq0;
    ray1.p = ray0.p - halfh * cnn0_csq0    * ray0.q;
    ray1.q = ray0.q + halfh * ccpx0.real() * ray0.p;
    
    // *** Phase 2
    
    EvaluateSSP(ray1.x, ccpx1, gradc1, crr1, crz1, czz1, rho, freq, ssp, iSegz, iSegr);
    csq1      = SQ(ccpx1.real());
    cnn1_csq1 = crr1 * SQ(ray1.t.y) - RC(2.0) * crz1 * ray1.t.x * ray1.t.y + czz1 * SQ(ray1.t.x);
    
    // The Munk test case with a horizontally launched ray caused problems.
    // The ray vertexes on an interface and can ping-pong around that interface.
    // Have to be careful in that case about big changes to the stepsize (that invalidate the leap-frog scheme) in phase II.
    // A modified Heun or Box method could also work.

    urayt1 = ccpx1.real() * ray1.t; // unit tangent
    
    // reduce h to land on boundary
    ReduceStep2D(ray0.x, urayt1, iSegz0, iSegr0, Topx, Topn, Botx, Botn, rTopSeg, rBotSeg,
        Beam, ssp, h, iSmallStepCtr);
    
    // use blend of f' based on proportion of a full step used.
    w1  = h / (RC(2.0) * halfh);
    w0  = RC(1.0) - w1;
    hw0 = h * w0;
    hw1 = h * w1;

    ray2->x   = ray0.x   + hw0 * urayt0                 + hw1 * urayt1;
    ray2->t   = ray0.t   - hw0 * gradc0 / csq0          - hw1 * gradc1 / csq1;
    ray2->p   = ray0.p   - hw0 * cnn0_csq0    * ray0.q - hw1 * cnn1_csq1    * ray1.q;
    ray2->q   = ray0.q   + hw0 * ccpx0.real() * ray0.p + hw1 * ccpx1.real() * ray1.p;
    ray2->tau = ray0.tau + hw0 / ccpx0                  + hw1 / ccpx1;

    ray2->Amp       = ray0.Amp;
    ray2->Phase     = ray0.Phase;
    ray2->NumTopBnc = ray0.NumTopBnc;
    ray2->NumBotBnc = ray0.NumBotBnc;
    
    // If we crossed an interface, apply jump condition

    EvaluateSSP(ray2->x, ccpx2, gradc2, crr2, crz2, czz2, rho, freq, ssp, iSegz, iSegr);
    ray2->c = ccpx2.real();
    
    if(iSegz != iSegz0 || iSegr != iSegr0){
        gradcjump = gradc2 - gradc0;
        ray2n = vec2(-ray2->t.y, ray2->t.x); // ray normal
        
        cnjump = glm::dot(gradcjump, ray2n);
        csjump = glm::dot(gradcjump, ray2->t);
        
        if(iSegz != iSegz0){             // crossing in depth
            rm =  ray2->t.x / ray2->t.y; // this is tan( alpha ) where alpha is the angle of incidence
        }else{                           // crossing in range
            rm = -ray2->t.y / ray2->t.x; // this is tan( alpha ) where alpha is the angle of incidence
        }                                // LP: The case where it crosses in depth and range simultaneously is not handled.
        
        rn = rm * (RC(2.0) * cnjump - rm * csjump) / ccpx2.real();
        ray2->p = ray2->p - ray2->q * rn;
    }
}
