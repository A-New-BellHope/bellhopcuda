#pragma once
#include "common.hpp"
#include "ldio.hpp"

inline void ReadPat(std::string FileRoot, std::ostream &PRTFile,
    int32_t &NSBPPts, real *&SrcBmPat, char SBPFlag)
{
    if(SBPFlag == '*'){
        PRTFile << "\n______________________________\nUsing source beam pattern file\n";
        
        LDIFile SBPFile(WithExtension(FileRoot, ".sbp"));
        if(!SBPFile.Good()){
            PRTFile << "SBPFile = " << WithExtension(FileRoot, ".sbp") << ".sbp\n";
            std::cout << "BELLHOP-ReadPat: Unable to open source beampattern file\n";
            std::abort();
        }
        
        SBPFile.List(); SBPFile.Read(NSBPPts);
        PRTFile << "Number of source beam pattern points " << NSBPPts << "\n";
        
        SrcBmPat = allocate<real>(NSBPPts * 2);
        
        PRTFile << "\n Angle (degrees)  Power (dB)\n" << std::setprecision(3);
        
        for(int32_t i=0; i<NSBPPts; ++i){
            SBPFile.List(); SBPFile.Read(&SrcBmPat[2*i], 2);
            PRTFile << SrcBmPat[2*i] << " " << SrcBmPat[2*i+1] << "\n";
        }
    }else{
        NSBPPts = 2;
        SrcBmPat = allocate<real>(2*2);
        SrcBmPat[0*2+0] = RC(-180.0); SrcBmPat[0*2+1] = RC(0.0);
        SrcBmPat[1*2+0] = RC( 180.0); SrcBmPat[1*2+1] = RC(0.0);
    }
    
    // convert dB to linear scale
    for(int32_t i=0; i<NSBPPts; ++i) SrcBmPat[i*2+1] = 
        STD::pow(RC(10.0), SrcBmPat[i*2+1] / RC(20.0));
}
