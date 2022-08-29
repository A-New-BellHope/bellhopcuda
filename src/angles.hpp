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

/**
 * LP: RunType was originally length 6 in the FORTRAN, but Beam->RunType was
 * being passed to it, which is definitely length 7.
 */
inline void ReadRayElevationAngles(real freq, real Depth,
    const char (&TopOpt)[6], const char (&RunType)[7],
    LDIFile &ENVFile, PrintFileEmu &PRTFile,
    AnglesStructure *Angles, Position *Pos)
{
    constexpr real c0 = FL(1500.0);
    
    if(TopOpt[5] == 'I'){
        // option to trace a single beam
        LIST(ENVFile); ENVFile.Read(Angles->Nalpha); ENVFile.Read(Angles->iSingle_alpha);
    }else{
        LIST(ENVFile); ENVFile.Read(Angles->Nalpha);
    }
    
    if(Angles->Nalpha == 0){ // automatically estimate Nalpha to use
        if(RunType[0] == 'R'){
            Angles->Nalpha = 50; // For a ray trace plot, we don't want too many rays ...
        }else{
            // you're letting ME choose? OK: ideas based on an isospeed ocean
            // limit based on phase of adjacent beams at maximum range
            Angles->Nalpha = bhc::max((int)(FL(0.3) * Pos->Rr[Pos->NRr-1] * freq / c0), 300);
            
            // limit based on having beams that are thin with respect to the water depth
            // assumes also a full 360 degree angular spread of rays
            // Should check which Depth is used here, in case where there is a variable bathymetry
            real d_theta_recommended = STD::atan(Depth / (FL(10.0) * Pos->Rr[Pos->NRr-1]));
            Angles->Nalpha = bhc::max((int)(REAL_PI / d_theta_recommended), Angles->Nalpha);
        }
    }
    
    checkallocate(Angles->alpha, bhc::max(3, Angles->Nalpha));
    
    if(Angles->Nalpha > 2) Angles->alpha[2] = FL(-999.9);
    LIST(ENVFile); ENVFile.Read(Angles->alpha, Angles->Nalpha);
    
    SubTab(Angles->alpha, Angles->Nalpha);
    Sort(  Angles->alpha, Angles->Nalpha);
    
    CheckFix360Sweep(Angles->alpha, Angles->Nalpha);
    
    PRTFile << "__________________________________________________________________________\n\n";
    PRTFile << "Number of beams in elevation   = " << Angles->Nalpha << "\n";
    if(Angles->iSingle_alpha > 0) PRTFile << "Trace only beam number " << Angles->iSingle_alpha << "\n";
    PRTFile << "Beam take-off angles (degrees)\n";
    
    EchoVector(Angles->alpha, Angles->Nalpha, PRTFile);
    
    if(Angles->Nalpha > 1 && Angles->alpha[Angles->Nalpha-1] == Angles->alpha[0]){
        GlobalLog("ReadRayElevationAngles: First and last beam take-off angle are identical\n");
        std::abort();
    }
    
    if(TopOpt[5] == 'I'){
        if(Angles->iSingle_alpha < 1 || Angles->iSingle_alpha > Angles->Nalpha){
            GlobalLog("ReadRayElevationAngles: Selected beam, iSingl not in [1, Angles->Nalpha]\n");
            std::abort();
        }
    }
}

}
