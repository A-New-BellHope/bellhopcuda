#pragma once
#include "common.hpp"
#include "ldio.hpp"

struct rxyz {
    real r, x, y, z;
};

struct BeamStructure {
    //LP: NSteps moved out of this struct as it's a property of a single beam.
    int32_t NBeams, Nimage, iBeamWindow;
    real deltas, epsMultiplier, rLoop;
    char Component;
    char Type[4];
    char RunType[7];
    rxyz Box;
};

struct BeamInfo {
    int32_t NSBPPts;
    real *SrcBmPat;
    char SBPFlag;
};

inline void ReadPat(std::string FileRoot, std::ofstream &PRTFile,
    BeamInfo *beaminfo)
{
    if(beaminfo->SBPFlag == '*'){
        PRTFile << "\n______________________________\nUsing source beam pattern file\n";
        
        LDIFile SBPFile(FileRoot, ".sbp");
        if(!SBPFile.Good()){
            PRTFile << "SBPFile = " << FileRoot << ".sbp\n";
            std::cout << "BELLHOP-ReadPat: Unable to open source beampattern file\n";
            std::abort();
        }
        
        SBPFile.List(); SBPFile.Read(beaminfo->NSBPPts);
        PRTFile << "Number of source beam pattern points " << beaminfo->NSBPPts << "\n";
        
        beaminfo->SrcBmPat = allocate<real>(beaminfo->NSBPPts * 2);
        
        PRTFile << "\n Angle (degrees)  Power (dB)\n" << std::setprecision(3);
        
        for(int32_t i=0; i<beaminfo->NSBPPts; ++i){
            SBPFile.List(); SBPFile.Read(&beaminfo->SrcBmPat[2*i], 2);
            PRTFile << beaminfo->SrcBmPat[2*i] << " " << beaminfo->SrcBmPat[2*i+1] << "\n";
        }
    }else{
        beaminfo->NSBPPts = 2;
        beaminfo->SrcBmPat = allocate<real>(2*2);
        beaminfo->SrcBmPat[0*2+0] = RC(-180.0); beaminfo->SrcBmPat[0*2+1] = RC(0.0);
        beaminfo->SrcBmPat[1*2+0] = RC( 180.0); beaminfo->SrcBmPat[1*2+1] = RC(0.0);
    }
    
    //LP: BUG: BELLHOP does not require that the angles are monotonically
    //increasing. However, in bellhop.f90:261, it looks for the largest angle
    //below the target angle and then linearly interpolates between that one and
    //the next angle in the array. This is only sensible if the angles are
    //monotonically increasing.
    if(!monotonic(beaminfo->SrcBmPat, beaminfo->NSBPPts, 2, 0)){
        std::cout << "Source beam pattern angles are not monotonic. This was not "
            << "a requirement in BELLHOP but has been added to bellhopcuda "
            << "because BELLHOP gave nonsense results if they were not.\n";
        std::abort();
    }
    
    // convert dB to linear scale
    for(int32_t i=0; i<beaminfo->NSBPPts; ++i) beaminfo->SrcBmPat[i*2+1] = 
        STD::pow(RC(10.0), beaminfo->SrcBmPat[i*2+1] / RC(20.0));
}

/**
 * Limits for tracing beams
 */
inline void ReadBeamInfo(LDIFile &ENVFile, std::ofstream &PRTFile,
    BeamStructure *Beam)
{
    PRTFile << "\n__________________________________________________________________________\n\n";

    ENVFile.List(); ENVFile.Read(Beam->deltas); ENVFile.Read(Beam->Box.z); ENVFile.Read(Beam->Box.r);
    
    PRTFile << std::setprecision(4);
    PRTFile << "\n Step length,       deltas = " << std::setw(11) << Beam->deltas << " m\n\n";
    PRTFile << "Maximum ray depth, Box.z  = " << std::setw(11) << Beam->Box.z << " m\n";
    PRTFile << "Maximum ray range, Box.r  = " << std::setw(11) << Beam->Box.r << "km\n";
    
    Beam->Box.r *= RC(1000.0); // convert km to m
    
    // *** Beam characteristics ***
    
    Beam->Type[3] = Beam->RunType[6]; // selects beam shift option
    
    if(Beam->Type[3] == 'S'){
        PRTFile << "Beam shift in effect\n";
    }else{
        PRTFile << "No beam shift in effect\n";
    }
    
    if(Beam->RunType[0] != 'R'){ // no worry about the beam type if this is a ray trace run
        
        // Beam%Type( 1 : 1 ) is
        //   'G" or "^' Geometric hat beams in Cartesian coordinates
        //   'g' Geometric hat beams in ray-centered coordinates
        //   'B' Geometric Gaussian beams in Cartesian coordinates
        //   'b' Geometric Gaussian beams in ray-centered coordinates
        //   'S' Simple Gaussian beams
        //   'C' Cerveny Gaussian beams in Cartesian coordinates
        //   'R' Cerveny Gaussian beams in Ray-centered coordinates
        // Beam%Type( 2 : 2 ) controls the setting of the beam width
        //   'F' space Filling
        //   'M' minimum width
        //   'W' WKB beams
        // Beam%Type( 3 : 3 ) controls curvature changes on boundary reflections
        //   'D' Double
        //   'S' Single
        //   'Z' Zero
        // Beam%Type( 4 : 4 ) selects whether beam shifts are implemented on boundary reflection
        //   'S' yes
        //   'N' no
 
        // Curvature change can cause overflow in grazing case
        // Suppress by setting BeamType( 3 : 3 ) = 'Z'
        
        Beam->Type[0] = Beam->RunType[1];
        switch(Beam->Type[0]){
        case 'G':
        case 'g':
        case '^':
        case 'B':
        case 'b':
        case 'S':
            break;
        case 'R':
        case 'C':
            ENVFile.List(); ENVFile.Read(&Beam->Type[1], 2); ENVFile.Read(Beam->epsMultiplier); ENVFile.Read(Beam->rLoop);
            PRTFile << "\n\nType of beam = " << Beam->Type[0] << "\n";
            switch(Beam->Type[2]){
            case 'D':
                PRTFile << "Curvature doubling invoked\n"; break;
            case 'Z':
                PRTFile << "Curvature zeroing invoked\n"; break;
            case 'S':
                PRTFile << "Standard curvature condition\n"; break;
            default:
                std::cout << "ReadEnvironment: Unknown curvature condition\n";
                std::abort();
            }
            break;
        default:
            std::cout << "ReadEnvironment: Unknown beam type (second letter of run type\n";
            std::abort();
        }
    }
}
