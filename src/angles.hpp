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
#include "runtype.hpp"

namespace bhc {

/**
 * LP: RunType now passed as part of Beam.
 */
template<bool O3D, bool BEARING> inline void ReadRayAngles(real freq, real Depth,
    const char (&TopOpt)[6], LDIFile &ENVFile, PrintFileEmu &PRTFile,
    AngleInfo &a, Position *Pos, const BeamStructure<O3D> *Beam)
{
    constexpr real c0 = FL(1500.0);
    const char *const FuncName = BEARING ? "ReadRayBearingAngles" : "ReadRayElevationAngles";
    
    if(TopOpt[5] == 'I'){
        // option to trace a single beam
        LIST(ENVFile); ENVFile.Read(a.n); ENVFile.Read(a.iSingle);
    }else{
        LIST(ENVFile); ENVFile.Read(a.n);
    }
    
    if(a.n == 0){ // automatically estimate n to use
        if(IsRayRun(Beam)){
            a.n = 50; // For a ray trace plot, we don't want too many rays ...
        }else{
            // you're letting ME choose? OK: ideas based on an isospeed ocean
            // limit based on phase of adjacent beams at maximum range
            a.n = bhc::max((int)((BEARING ? FL(0.1) : FL(0.3)) * Pos->Rr[Pos->NRr-1] * freq / c0), 300);
            
            if constexpr(!BEARING){
                // limit based on having beams that are thin with respect to the water depth
                // assumes also a full 360 degree angular spread of rays
                // Should check which Depth is used here, in case where there is a variable bathymetry
                real d_theta_recommended = STD::atan(Depth / (FL(10.0) * Pos->Rr[Pos->NRr-1]));
                a.n = bhc::max((int)(REAL_PI / d_theta_recommended), a.n);
            }
        }
    }
    
    checkallocate(a.angles, bhc::max(3, a.n));
    
    if(a.n > 2) a.angles[2] = FL(-999.9);
    LIST(ENVFile); ENVFile.Read(a.angles, a.n);
    
    SubTab(a.angles, a.n);
    Sort(  a.angles, a.n);
    
    CheckFix360Sweep(a.angles, a.n);
    
    if constexpr(BEARING){
        // Nx2D CASE: beams must lie on rcvr radials--- replace beta with theta
        if(Beam->RunType[5] == '2' && !IsRayRun(Beam)){
            PRTFile << "\nReplacing beam take-off angles, beta, with receiver bearing lines, theta\n";
            checkdeallocate(a.angles);
            
            a.n = Pos->Ntheta;
            checkallocate(a.angles, bhc::max(3, a.n));
            for(int32_t i=0; i<a.n; ++i) a.angles[i] = Pos->theta[i]; // a.n should = Pos->Ntheta
        }
    }
    
    if constexpr(!BEARING){
        PRTFile << "__________________________________________________________________________\n";
    }
    PRTFile << "\n   Number of beams in " << (BEARING ? "bearing  " : "elevation") << "   = " << a.n << "\n";
    if(a.iSingle > 0) PRTFile << "Trace only beam number " << a.iSingle << "\n";
    PRTFile << "   Beam take-off angles (degrees)\n";
    
    EchoVector(a.angles, a.n, PRTFile);
    
    if(a.n > 1 && a.angles[a.n-1] == a.angles[0]){
        GlobalLog("%s: First and last beam take-off angle are identical\n", FuncName);
        std::abort();
    }
    
    if(TopOpt[5] == 'I'){
        if(a.iSingle < 1 || a.iSingle > a.n){
            GlobalLog("%s: Selected beam, iSingle not in [1, a.n]\n", FuncName);
            std::abort();
        }
    }
    
    // LP: This and the d computation below was in setup / BellhopCore for alpha,
    // but here for beta. 
    for(int32_t i=0; i<a.n; ++i)
        a.angles[i] *= DegRad; // convert to radians
    
    a.d = FL(0.0);
    if(a.n != 1)
        a.d = (a.angles[a.n-1] - a.angles[0]) / (a.n-1); // angular spacing between beams
    
}

}
