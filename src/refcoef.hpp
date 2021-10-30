#pragma once
#include "common.hpp"
#include "ldio.hpp"

struct ReflectionCoef {
    real theta, r, phi;
};

struct ReflectionInfo {
    int32_t NBotPts, NTopPts;
    ReflectionCoef *RBot, *RTop;
};

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
    real alpha, thetaIntr;
    
    iLeft = 0;
    iRight = NPts - 1;
    
    thetaIntr = RInt.theta.real(); // This should be unnecessary? probably used when I was doing complex angles
    
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

/**
 * Optionally read in reflection coefficient for Top or Bottom boundary
 * 
 * BotRC, TopRC: flag set to 'F' if refl. coef. is to be read from a File
 */
inline void ReadReflectionCoefficient(std::string FileRoot, char BotRC, char TopRC,
    std::ofstream &PRTFile, ReflectionInfo *refl)
{
    if(BotRC == 'F'){
        PRTFile << "__________________________________________________________________________\n\n";
        PRTFile << "Using tabulated bottom reflection coef.\n";
        LDIFile BRCFile(FileRoot + ".brc");
        if(!BRCFile.Good()){
            PRTFile << "BRCFile = " << FileRoot + ".brc\n";
            std::cout << "ReadReflectionCoefficient: Unable to open Bottom Reflection Coefficient file\n";
            std::abort();
        }
        
        BRCFile.List(); BRCFile.Read(refl->NBotPts);
        PRTFile << "Number of points in bottom reflection coefficient = " << refl->NBotPts << "\n";
        
        if(refl->RBot != nullptr) deallocate(refl->RBot);
        refl->RBot = allocate<ReflectionCoef>(refl->NBotPts);
        
        BRCFile.List();
        for(int32_t itheta=0; itheta<refl->NBotPts; ++itheta){
            BRCFile.Read(refl->RBot[itheta].theta);
            BRCFile.Read(refl->RBot[itheta].r);
            BRCFile.Read(refl->RBot[itheta].phi);
            refl->RBot->phi *= DegRad; // convert to radians
        }
    }else{ // should allocate something anyway, since variable is passed
        //LP: BUG: BELLHOP does not deallocate this array first, despite doing
        //so above.
        if(refl->RBot != nullptr) deallocate(refl->RBot);
        refl->RBot = allocate<ReflectionCoef>(1);
    }
    
    // Optionally read in top reflection coefficient
    
    if(TopRC == 'F'){
        PRTFile << "__________________________________________________________________________\n\n";
        PRTFile << "Using tabulated top    reflection coef.\n";
        LDIFile TRCFile(FileRoot + ".trc");
        if(!TRCFile.Good()){
            PRTFile << "TRCFile = " << FileRoot + ".trc\n";
            std::cout << "ReadReflectionCoefficient: Unable to open Top Reflection Coefficient file\n";
            std::abort();
        }
        
        TRCFile.List(); TRCFile.Read(refl->NTopPts);
        PRTFile << "Number of points in top reflection coefficient = " << refl->NTopPts << "\n";
        
        if(refl->RTop != nullptr) deallocate(refl->RTop);
        refl->RTop = allocate<ReflectionCoef>(refl->NTopPts);
        
        TRCFile.List();
        for(int32_t itheta=0; itheta<refl->NTopPts; ++itheta){
            TRCFile.Read(refl->RTop[itheta].theta);
            TRCFile.Read(refl->RTop[itheta].r);
            TRCFile.Read(refl->RTop[itheta].phi);
            refl->RTop->phi *= DegRad; // convert to radians
        }
    }else{ // should allocate something anyway, since variable is passed
        //LP: BUG: BELLHOP does not deallocate this array first, despite doing
        //so above.
        if(refl->RTop != nullptr) deallocate(refl->RTop);
        refl->RTop = allocate<ReflectionCoef>(1);
    }
    
    // Optionally read in internal reflection coefficient data
    
    if(BotRC == 'P'){
        std::cout << "Internal reflections not supported by BELLHOP and therefore "
            "not supported by bellhopcuda\n";
        std::abort();
    }
}
