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

namespace bhc {

inline void ReadPat(std::string FileRoot, PrintFileEmu &PRTFile,
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
        
        checkallocate(beaminfo->SrcBmPat, beaminfo->NSBPPts * 2);
        
        PRTFile << "\n Angle (degrees)  Power (dB)\n" << std::setprecision(3);
        
        for(int32_t i=0; i<beaminfo->NSBPPts; ++i){
            LIST(SBPFile); SBPFile.Read(&beaminfo->SrcBmPat[2*i], 2);
            PRTFile << beaminfo->SrcBmPat[2*i] << " " << beaminfo->SrcBmPat[2*i+1] << "\n";
        }
    }else{
        beaminfo->NSBPPts = 2;
        checkallocate(beaminfo->SrcBmPat, 2*2);
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
inline void ReadBeamInfo(LDIFile &ENVFile, PrintFileEmu &PRTFile,
    BeamStructure *Beam, const BdryType *Bdry, bool o3d)
{
    if(o3d){
        LIST(ENVFile); ENVFile.Read(Beam->deltas); ENVFile.Read(Beam->Box.x);
        ENVFile.Read(Beam->Box.y); ENVFile.Read(Beam->Box.z);
        Beam->Box.x *= FL(1000.0); // convert km to m
        Beam->Box.y *= FL(1000.0); // convert km to m
        
        if(Beam->deltas == FL(0.0)) Beam->deltas = (Bdry->Bot.hs.Depth - Bdry->Top.hs.Depth) / FL(10.0); // Automatic step size selection
    }else{
        LIST(ENVFile); ENVFile.Read(Beam->deltas); ENVFile.Read(Beam->Box.z); ENVFile.Read(Beam->Box.r);
    }
    
    PRTFile << std::setprecision(4);
    PRTFile << "\n Step length,       deltas = " << std::setw(11) << Beam->deltas << " m\n\n";
    if(o3d){
        PRTFile << "Maximum ray x-range, Box.x  = " << std::setw(11) << Beam->Box.x << " m\n";
        PRTFile << "Maximum ray y-range, Box.x  = " << std::setw(11) << Beam->Box.y << " m\n";
        PRTFile << "Maximum ray z-range, Box.x  = " << std::setw(11) << Beam->Box.z << " m\n";
    }else{
        PRTFile << "Maximum ray depth, Box.z  = " << std::setw(11) << Beam->Box.z << " m\n";
        PRTFile << "Maximum ray range, Box.r  = " << std::setw(11) << Beam->Box.r << "km\n";
        
        Beam->Box.r *= FL(1000.0); // convert km to m
    }
    
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

}
