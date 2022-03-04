#pragma once
#include "common.hpp"

inline void ReadPat(std::string FileRoot, std::ostream &PRTFile,
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
        
        LIST(SBPFile); SBPFile.Read(beaminfo->NSBPPts);
        PRTFile << "Number of source beam pattern points " << beaminfo->NSBPPts << "\n";
        
        beaminfo->SrcBmPat = allocate<real>(beaminfo->NSBPPts * 2);
        
        PRTFile << "\n Angle (degrees)  Power (dB)\n" << std::setprecision(3);
        
        for(int32_t i=0; i<beaminfo->NSBPPts; ++i){
            LIST(SBPFile); SBPFile.Read(&beaminfo->SrcBmPat[2*i], 2);
            PRTFile << beaminfo->SrcBmPat[2*i] << " " << beaminfo->SrcBmPat[2*i+1] << "\n";
        }
    }else{
        beaminfo->NSBPPts = 2;
        beaminfo->SrcBmPat = allocate<real>(2*2);
        beaminfo->SrcBmPat[0*2+0] = FL(-180.0); beaminfo->SrcBmPat[0*2+1] = FL(0.0);
        beaminfo->SrcBmPat[1*2+0] = FL( 180.0); beaminfo->SrcBmPat[1*2+1] = FL(0.0);
    }
    
    if(!monotonic(beaminfo->SrcBmPat, beaminfo->NSBPPts, 2, 0)){
        std::cout << "BELLHOP-ReadPat: Source beam pattern angles are not monotonic\n";
        std::abort();
    }
    
    // convert dB to linear scale
    for(int32_t i=0; i<beaminfo->NSBPPts; ++i) beaminfo->SrcBmPat[i*2+1] = 
        STD::pow(FL(10.0), beaminfo->SrcBmPat[i*2+1] / FL(20.0));
}

/**
 * Limits for tracing beams
 */
inline void ReadBeamInfo(LDIFile &ENVFile, std::ostream &PRTFile,
    BeamStructure *Beam)
{
    PRTFile << "\n__________________________________________________________________________\n\n";

    LIST(ENVFile); ENVFile.Read(Beam->deltas); ENVFile.Read(Beam->Box.z); ENVFile.Read(Beam->Box.r);
    
    PRTFile << std::setprecision(4);
    PRTFile << "\n Step length,       deltas = " << std::setw(11) << Beam->deltas << " m\n\n";
    PRTFile << "Maximum ray depth, Box.z  = " << std::setw(11) << Beam->Box.z << " m\n";
    PRTFile << "Maximum ray range, Box.r  = " << std::setw(11) << Beam->Box.r << "km\n";
    
    Beam->Box.r *= FL(1000.0); // convert km to m
    
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
            LIST(ENVFile); ENVFile.Read(&Beam->Type[1], 2); ENVFile.Read(Beam->epsMultiplier); ENVFile.Read(Beam->rLoop);
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
            
            PRTFile << "Epsilon multiplier " << Beam->epsMultiplier << "\n";
            PRTFile << "Range for choosing beam width " << Beam->rLoop << "\n";
            
            // Images, windows
            // LP: These values are not initialized if not written in the file,
            // and Component is not always written in the test env files.
            Beam->Nimage = 1;
            Beam->iBeamWindow = 4;
            Beam->Component = 'P';
            LIST(ENVFile); ENVFile.Read(Beam->Nimage); ENVFile.Read(Beam->iBeamWindow); ENVFile.Read(Beam->Component);
            PRTFile << "\nNumber of images, Nimage  = " << Beam->Nimage << "\n";
            PRTFile << "Beam windowing parameter  = " << Beam->iBeamWindow << "\n";
            PRTFile << "Component                 = " << Beam->Component << "\n";
            break;
        default:
            std::cout << "ReadEnvironment: Unknown beam type (second letter of run type\n";
            std::abort();
        }
        
        PRTFile << "\n";
    }
}
