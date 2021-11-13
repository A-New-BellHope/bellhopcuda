#pragma once
#include "common.hpp"

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
