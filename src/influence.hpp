#pragma once
#include "common.hpp"
#include "atomics.hpp"

struct InfluenceRayInfo {
    real phase, qOld, q0;
};

/**
 * Calculates a smoothing function based on the h0 hermite cubic
 * x is the point where the function is to be evaluated
 * returns:
 * [ 0, x1 ] = 1
 * [x1, x2 ] = cubic taper from 1 to 0
 * [x2, inf] = 0
 */
HOST_DEVICE inline real Hermite(real x, real x1, real x2)
{
    real Ax = STD::abs(x);
    if(Ax <= x1){
        return RC(1.0);
    }else if(Ax >= x2){
        return RC(0.0);
    }else{
        real u = (Ax - x1) / (x2 - x1);
        return (RC(1.0) + RC(2.0) * u) * SQ(RC(1.0) - u);
    }
    // ret /= (RC(0.5) * (x1 + x2));
}

/**
 * Scale the pressure field
 * 
 * r: ranges (LP: [Nr])
 * Dalpha: angular spacing between rays
 * freq: source frequency
 * c: nominal sound speed (LP: [NRz][Nr])
 * u: Pressure field
 */
HOST_DEVICE inline void ScalePressure(real Dalpha, real c, real *r, 
    cpx *u, int32_t NRz, int32_t Nr, char (&RunType)[5], real freq)
{
    // Compute scale factor for field
    real cnst;
    if(RunType[1] == 'C' || RunType[1] == 'R'){
        // Cerveny Gaussian beams in Cartesian or Ray-centered coordinates
        cnst = -Dalpha * STD::sqrt(freq) / c;
    }else{
        cnst = RC(-1.0);
    }
    
    // For incoherent run, convert intensity to pressure
    if(RunType[0] != 'C'){
        for(int32_t irz=0; irz<NRz; ++irz){
            for(int32_t ir=0; ir<Nr; ++ir){
                u[irz*Nr+ir] = cpx(STD::sqrt(u[irz*Nr+ir].real()), RC(0.0));
            }
        }
    }
    
    // scale and/or incorporate cylindrical spreading
    for(int32_t ir=0; ir<Nr; ++ir){
        real factor;
        if(RunType[3] == 'X'){ // line source
            factor = RC(-4.0) * STD::sqrt(M_PI) * cnst;
        }else{ // point source
            if(r[ir] == RC(0.0)){
                factor = RC(0.0); // avoid /0 at origin, return pressure = 0
            }else{
                factor = cnst / STD::sqrt(STD::abs(r[ir]));
            }
        }
        for(int32_t irz=0; irz<NRz; ++irz){
            u[irz*Nr+ir] *= factor;
        }
    }
}

/**
 * Checks for a branch cut crossing and updates kmah accordingly
 */
HOST_DEVICE inline void BranchCut(const cpx &q1C, const cpx &q2C,
    char (&BeamType)[4], int32_t &kmah)
{
    real q1, q2;
    if(BeamType[1] == 'W'){ // WKBeams
        q1 = q1C.real();
        q2 = q2C.real();
    }else{
        if(q2C.real() >= RC(0.0)) return;
        q1 = q1C.imag();
        q2 = q2C.imag();
    }
    if( (q1 < RC(0.0) && q2 >= RC(0.0)) || 
        (q1 > RC(0.0) && q2 <= RC(0.0)) ) kmah = -kmah;
}

/**
 * phase shifts at caustics
 */
HOST_DEVICE inline void IncPhaseIfCaustic(InfluenceRayInfo &inflray, real q)
{
    if(q < RC(0.0) && inflray.qOld >= RC(0.0) || q > RC(0.0) && inflray.qOld <= RC(0.0))
        inflray.phase += M_PI * RC(0.5);
}

HOST_DEVICE inline void Init_InfluenceSGB(InfluenceRayInfo &inflray)
{
    inflray.phase = RC(0.0);
    inflray.qOld = RC(1.0);
}

/**
 * Bucker's Simple Gaussian Beams in Cartesian coordinates
 */
HOST_DEVICE inline void Step_InfluenceSGB(
    const ray2DPt &point0, const ray2DPt &point1, InfluenceRayInfo &inflray,
    int32_t is,
    cpx *u, real alpha, real dalpha, real RadiusMax,
    int32_t NRz_per_range, const Position *Pos, const BeamStructure *Beam)
{
    real w;
    vec2 x, rayt;
    cpx tau;
    
    real Ratio1 = STD::sqrt(STD::cos(alpha));
    const real beta = RC(0.98); // Beam Factor
    real a = RC(-4.0) * STD::log(beta) / SQ(dalpha);
    real cn = Dalpha * STD::sqrt(a / M_PI);
    real rA = point0.x.x;
    int32_t ir = 0;
    
    real rB = point1.x.x;
    
    real q = point0.q.x;
    IncPhaseIfCaustic(inflray, q);
    inflray.qOld = q;
    
    // Loop over bracketted receiver ranges
    // LP: BUG: This does not check for "bracketted" receivers, only that the
    // end of the ray segment is to the right of the receiver. w may be an
    // arbitrarily low negative number, which can lead to garbage extrapolated
    // values for the other variables. See the BUG below.
    while(STD::abs(rB - rA) > RC(1.0e3) * spacing(rA) && rB > Pos->Rr[ir]){
        w    = (Pos->Rr[ir] - rA) / (rB - rA);
        x    = point0.x   + w * (point1.x   - point0.x);
        rayt = point0.t   + w * (point1.t   - point0.t);
        q    = point0.q.x + w * (point1.q.x - point0.q.x);
        tau  = point0.tau + w * (point1.tau - point0.tau);
        
        // following is incorrect because ray doesn't always use a step of deltas
        // LP: The while ignores extremely small steps, but those small steps
        // still increment is, so the later ray segments still treat it as if
        // all steps leading up to them were of size deltas.
        real sint = ((real)is + w) * Beam->deltas;
        
        IncPhaseIfCaustic(inflray, q);
        
        for(int32_t iz=0; i<NRz_per_range; ++i){
            real deltaz = Pos->Rz[iz] - x.y; // ray to rcvr distance
            //LP: This is commented out, but seems very important to have.
            //Though with all the other bugs, perhaps this is not the main issue.
            //real Adeltaz = STD::abs(deltaz);
            //if(Adeltaz < RadiusMax){
                if(Beam->RunType[0] == 'E'){ // eigenrays
                    print("Eigenrays not yet supported\n");
                    bail();
                    //real SrcDeclAngle = RadDeg * alpha; // take-off angle in degrees
                    //WriteRay2D(SrcDeclAngle, is);
                }else{ // coherent TL
                    // LP: BUG: It may be incoherent or semi-coherent.
                    real cpa = STD::abs(deltaz * (rB - rA)) / STD::sqrt(SQ(rB - rQ) 
                        + SQ(point1.x.y - point0.x.y));
                    real ds = STD::sqrt(SQ(deltaz) - SQ(cpa));
                    real sx1 = sint + ds;
                    real thet = STD::atan(cpa / sx1);
                    real delay = tau + rayt.y * deltaz;
                    cpx contri = Ratio1 * cn * point1.Amp * STD::exp(-a * SQ(thet) -
                        J * (omega * delay - point1.Phase - inflray.phase)) / STD::sqrt(sx1);
                    AtomicAddCpx(&u[iz*Pos->NRr + ir], contri);
                }
            //}
        }
        
        inflray.qOld = q;
        ++ir;
        // LP: BUG: In the FORTRAN this is also return, but it is within the
        // loop over the ray steps. This means that once the ray is to the right
        // of all receivers, it ignores the rest of the ray. This is wrong for
        // two reasons. First, it assumes rays always travel to the right, which
        // is not always true for certain bathymetries. Second, the influence is
        // only computed when the ray is to the right of a receiver. If there
        // are receviers at different range positions, the left ones will
        // receive contributions for every step the ray is to the right. How
        // many of these steps happen depends on where the right receivers are,
        // which means the computed signal at one receiver depends on the
        // position of another receiver, which is absurd.
        if(ir >= Pos->NRr) return;
    }
}

HOST_DEVICE inline void ApplyContribution(real Amp, real cnst, real w,
    real omega, real delay, real phaseInt, cpx delay, 
    cpx *u, const BeamStructure *Beam)
{
    switch(Beam->RunType[0]){
    case 'E':
        // eigenrays
        print("Eigenrays not yet supported\n");
        bail();
        //WriteRay2D(SrcDeclAngle, is);
        break;
    case 'A':
    case 'a':
        // arrivals
        print("Arrivals not yet supported\n");
        bail();
        //AddArr(omega, iz, ir, Amp, phaseInt, delay, SrcDeclAngle, RcvrDeclAngle,
        //    point1.NumTopBnc, point1.NumBotBnc);
        break;
    case 'C':
        // coherent TL
        AtomicAddCpx(u, Amp * STD::exp(-J * (omega * delay - phaseInt)));
        // omega * SQ(n) / (RC(2.0) * SQ(point1.c) * delay)))) // curvature correction
        break;
    default:
        // incoherent/semicoherent TL
        real v = STD::pow(cnst * STD::exp((omega * delay).imag()), RC(2.0) * w);
        if(Beam->Type[0] == 'B'){
            // Gaussian beam
            v *= STD::sqrt(RC(2.0) * M_PI);
        }
        AtomicAddCpx(u, cpx(v, RC(0.0)));
    }
}

HOST_DEVICE inline void Init_InfluenceGeoHatOrGaussianCart(
    InfluenceRayInfo &inflray, const ray2DPt &point0, real dalpha)
{
    inflray.q0 = point0.c / Dalpha; // Reference for J = q0 / q
    inflray.phase = RC(0.0);
    inflray.qOld = point0.q.x; // used to track KMAH index
}

/**
 * Geometric, Gaussian beams in Cartesian coordintes
 *
 * alpha: take-off angle
 * dalpha: angular spacing
 * u: complex pressure field
 */
HOST_DEVICE inline void Step_InfluenceGeoHatOrGaussianCart(
    bool isGaussian, const ray2DPt &point0, const ray2DPt &point1,
    InfluenceRayInfo &inflray, int32_t is, cpx *u, real alpha, real dalpha,
    int32_t NRz_per_range, const Position *Pos, const BeamStructure *Beam)
{
    const int32_t BeamWindow = 4; // beam window: kills beams outside e**(-0.5 * ibwin**2 )
    
    vec2 x_ray, rayt, rayn, x_rcvr;
    real q0, SrcDeclAngle, RcvrDeclAngle, Ratio1;
    real rlen, RadiusMax, zMin, zMax, sigma, lambda, a, dqds;
    cpx dtauds;
    
    SrcDeclAngle = RadDeg * alpha;
    rA = point0.x.x;
    rB = point1.x.x;
    
    // what if never satistified?
    // what if there is a single receiver (ir = -1 possible)
    // LP: This implementation has been adjusted to handle these cases.
    int32_t ir = BinarySearchGEQ(Pos->Rr, Pos->NRr, 1, 0, STD::min(rA, rB));
    if(Pos->Rr[ir] < STD::min(rA, rB)) return; // not "bracketted"
    
    if(Beam->RunType[3] == 'R'){
        Ratio1 = STD::sqrt(STD::abs(STD::cos(alpha))); // point source
    }else{
        Ratio1 = RC(1.0); // line source
    }
    // sqrt( 2 * pi ) represents a sum of Gaussians in free space
    if(isGaussian){
        Ratio1 /= STD::sqrt(RC(2.0) * M_PI);
    }
    
    x_ray = point0.x;
    
    // compute normalized tangent (compute it because we need to measure the step length)
    rayt = point1.x - point0.x;
    rlen = glm::norm(tayt);
    if(rlen < RC(1.0e3) * spacing(point1.x.x)) return; // if duplicate point in ray, skip to next step along the ray
    rayt /= rlen;
    rayn = vec2(-rayt.y, rayt.x); // unit normal to ray
    RcvrDeclAngle = RadDeg * STD::atan2(rayt.y, rayt.x);
    
    dqds   = point1.q.x - point0.q.x;
    dtauds = point1.tau - point0.tau;
    
    q = point0.q.x;
    IncPhaseIfCaustic(inflray, q);
    inflray.qOld = q;
    
    sigma = STD::max(STD::abs(point0.q.x), STD::abs(point1.q.x) / 
        (inflray.q0 * STD::abs(rayt.x)); // beam radius projected onto vertical line
    if(isGaussian){
        // calculate beam width
        lambda    = point0.c / freq;
        sigma     = STD::max(sigma, STD::min(RC(0.2) * freq * point1.tau, M_PI * lambda));
        RadiusMax = BeamWindow * sigma;
    }else{
        RadiusMax = sigma;
    }
    
    // depth limits of beam
    if(STD::abs(rayt.x) > RC(0.5)){ // shallow angle ray
        zmin = min(point0.x.y, point1.x.y) - RadiusMax;
        zmax = max(point0.x.y, point1.x.y) + RadiusMax;
    }else{
        zmin = -REAL_MAX;
        zmax =  REAL_MAX;
    }
    
    // compute beam influence for this segment of the ray
    for(; ir<Pos->NRr; ++ir){
        // is Rr[ir] contained in [rA, rB)? Then compute beam influence
        // LP: Because of the new setup and always incrementing regardless of
        // which direction the ray goes, we only have to check this side.
        if(Pos->Rr[ir] >= STD::max(rA, rB)) return;
        
        for(int32_t iz=0; iz<NRz_per_range; ++iz){
            if(Beam->RunType[4] == 'I'){
                x_rcvr = vec2(Pos->Rr[ir], Pos->Rz[ir]); // irregular grid
            }else{
                x_rcvr = vec2(Pos->Rr[ir], Pos->Rz[iz]); // rectilinear grid
            }
            if(x_rcvr.y < zmin || x_rcvr.x > zmax) continue;
            
            s     = glm::dot(x_rcvr - x_ray, rayt) / rlen; // proportional distance along ray
            n     = STD::abs(glm::dot(x_rcvr - x_ray, rayn)); // normal distance to ray
            q     = point0.q.x + s * dqds; // interpolated amplitude
            sigma = STD::abs(q / inflray.q0);
            real beamWCompare;
            if(isGaussian){
                sigma = STD::max(sigma, STD::min(RC(0.2) * freq * point1.tau, M_PI * lambda)); // min pi * lambda, unless near
                beamWCompare = BeamWindow * sigma;
            }else{
                RadiusMax = sigma;
                beamWCompare = RadiusMax;
            }
            
            if(n < beamWCompare){ // Within beam window?
                delay = point0.tau + s * dtauds; // interpolated delay
                cnst  = Ratio1 * STD::sqrt(point1.c / STD::abs(q)) * point1.Amp;
                if(isGaussian){
                    a = STD::abs(inflray.q0 / q);
                    w = STD::exp(RC(-0.5) * SQ(n / sigma)) / (sigma * a); // Gaussian decay
                }else{
                    w = (RadiusMax - n) / RadiusMax; // hat function: 1 on center, 0 on edge
                }
                Amp      = cnst * w;
                phaseInt = (isGaussian ? point1.Phase : point0.Phase) + inflray.phase;
                // LP: BUG: point0/1.phase is discarded if this condition is met.
                if(q <= RC(0.0) && inflray.qOld > RC(0.0) || q >= RC(0.0) && inflray.qOld < RC(0.0))
                    phaseInt = inflray.phase + M_PI * RC(0.5); // phase shifts at caustics
                
                ApplyContribution(Amp, cnst, w, RC(2.0) * M_PI * freq0, 
                    delay, phaseInt, delay, &u[iz*Pos->NRr+ir], Beam);
            }
        }
    }
}
