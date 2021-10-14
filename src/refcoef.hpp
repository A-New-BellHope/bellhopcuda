#pragma once
#include "structs.hpp"

/**
 * Given an angle RInt%ThetaInt, returns the magnitude and
 * phase of the reflection coefficient (RInt%R, RInt%phi).
 * 
 * Uses linear interpolation using the two nearest abscissas
 * Assumes phi has been unwrapped so that it varies smoothly.
 * I tried modifying it to allow a complex angle of incidence but
 * stopped when I realized I needed to fuss with a complex atan2 routine
 * LP: C++ supports complex atan, though there is no complex atan2.
 * 
 * RInt: interpolated value of refl. coef.
 * r: Reflection coefficient table
 * NPts: # pts in refl. coef.
 */
HOST_DEVICE inline void InterpolateReflectionCoefficient(ReflectionCoef &RInt,
    const ReflectionCoef *r, int32_t NPts)
{
    int32_t iLeft, iRight, iMid;
    real alpha, Thetaintr;
    
    iLeft = 0;
    iRight = NPts - 1;
    
    thetaIntr = RInt.Theta.real(); // This should be unnecessary? probably used when I was doing complex angles
    
    // Three cases: ThetaInt left, in, or right of tabulated interval
    
    if(thetaIntr < r[iLeft].theta){
        // iRight = 1;
        RInt.r   = RC(0.0); // r[iLeft].r
        RInt.phi = RC(0.0); // r[iLeft].phi
        printf("Warning in InterpolateReflectionCoefficient : Refl. Coef. being "
            "set to 0 outside tabulated domain : angle = %f, lower limit = %f",
            thetaIntr, r[iLeft].theta);
    }else if(thetaIntr > r[iRight].theta){
        // iLeft = NPts - 2;
        RInt.r   = RC(0.0); // r[iRight].r
        RInt.phi = RC(0.0); // r[iRight].phi
        // printf("Warning in InterpolateReflectionCoefficient : Refl. Coef. being "
        //     "set to 0 outside tabulated domain : angle = %f, lower limit = %f",
        //     thetaIntr, r[iRight].theta);
    }else{
        // Search for bracketting abscissas: Log2( NPts ) stabs required for a bracket
        
        while(iLeft != iRight - 1){
            iMid = (iLeft + iRight) / 2;
            if(r[iMid].theta > thetaIntr){
                iRight = iMid;
            }else{
                iLeft = iMid;
            }
        }
        
        // Linear interpolation for reflection coef
        
        alpha    = (RInt.theta - r[iLeft].theta) / (r[iRight].theta - r[iLeft].theta);
        RInt.r   = (RC(1.0) - alpha) * r[iLeft].r   + alpha * r[iRight].r;
        RInt.phi = (RC(1.0) - alpha) * r[iLeft].phi + alpha * r[iRight].phi;
    }
}
