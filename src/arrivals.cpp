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
#include "arrivals.hpp"

namespace bhc {

void FinalizeArrivalsMode(const ArrInfo *arrinfo, const Position *Pos,
    const FreqInfo *freqinfo, const BeamStructure *Beam, std::string FileRoot, bool ThreeD)
{
    bool isAscii;
    LDOFile AARRFile;
    UnformattedOFile BARRFile;
    switch(Beam->RunType[0]){
    case 'A': // arrivals calculation, ascii
        isAscii = true;
        
        AARRFile.open(FileRoot + ".arr");
        AARRFile << (ThreeD ? "3D" : "2D") << '\n';
        AARRFile << freqinfo->freq0 << '\n';
        
        // write source locations
        if(ThreeD){
            AARRFile << Pos->NSx; AARRFile.write(Pos->Sx, Pos->NSx); AARRFile << '\n';
            AARRFile << Pos->NSy; AARRFile.write(Pos->Sy, Pos->NSy); AARRFile << '\n';
        }
        AARRFile << Pos->NSz; AARRFile.write(Pos->Sz, Pos->NSz); AARRFile << '\n';
        
        // write receiver locations
        AARRFile << Pos->NRz; AARRFile.write(Pos->Rz, Pos->NRz); AARRFile << '\n';
        AARRFile << Pos->NRr; AARRFile.write(Pos->Rr, Pos->NRr); AARRFile << '\n';
        if(ThreeD){
            AARRFile << Pos->Ntheta; AARRFile.write(Pos->theta, Pos->Ntheta); AARRFile << '\n';
        }
        break;
    case 'a': // arrivals calculation, binary
        isAscii = false;
        
        BARRFile.open(FileRoot + ".arr");
        BARRFile.rec(); BARRFile.write((ThreeD ? "'3D'" : "'2D'"), 4);
        BARRFile.rec(); BARRFile.write((float)freqinfo->freq0);
        
        // write source locations
        if(ThreeD){
            BARRFile.rec(); BARRFile.write(Pos->NSx); BARRFile.write(Pos->Sx, Pos->NSx);
            BARRFile.rec(); BARRFile.write(Pos->NSy); BARRFile.write(Pos->Sy, Pos->NSy);
        }
        BARRFile.rec(); BARRFile.write(Pos->NSz); BARRFile.write(Pos->Sz, Pos->NSz);
        
        // write receiver locations
        BARRFile.rec(); BARRFile.write(Pos->NRz); BARRFile.write(Pos->Rz, Pos->NRz);
        BARRFile.rec(); BARRFile.write(Pos->NRr); BARRFile.write(Pos->Rr, Pos->NRr);
        if(ThreeD){
            BARRFile.rec(); BARRFile.write(Pos->Ntheta); BARRFile.write(Pos->theta, Pos->Ntheta);
        }
        break;
    default:
        GlobalLog("FinalizeArrivalsMode called while not in arrivals mode\n");
        bail();
        return;
    }
    for(int32_t isrc = 0; isrc < Pos->NSz; ++isrc){
        int32_t maxn = 0;
        for(int32_t iz=0; iz<Pos->NRz_per_range; ++iz){
            for(int32_t ir=0; ir<Pos->NRr; ++ir){
                int32_t base = (isrc * Pos->NRz_per_range + iz) * Pos->NRr + ir;
                if(arrinfo->NArr[base] > maxn) maxn = arrinfo->NArr[base];
            }
        }
        if(isAscii){
            AARRFile << maxn << '\n';
        }else{
            BARRFile.rec(); BARRFile.write(maxn);
        }
        
        for(int32_t iz=0; iz<Pos->NRz_per_range; ++iz){
            for(int32_t ir=0; ir<Pos->NRr; ++ir){
                int32_t base = (isrc * Pos->NRz_per_range + iz) * Pos->NRr + ir;
                real factor;
                if(Beam->RunType[3] == 'X'){ // line source
                    factor = FL(4.0) * STD::sqrt(REAL_PI);
                }else{                       // point source
                    if(Pos->Rr[ir] == FL(0.0)){
                        factor = FL(1e5); // avoid /0 at origin
                    }else{
                        factor = FL(1.0) / STD::sqrt(Pos->Rr[ir]); // cyl. spreading
                    }
                }
                
                int32_t narr = arrinfo->NArr[base];
                if(isAscii){
                    AARRFile << narr << '\n';
                }else{
                    BARRFile.rec(); BARRFile.write(narr);
                }
                
                for(int32_t iArr=0; iArr<narr; ++iArr){
                    Arrival *arr = &arrinfo->Arr[base * arrinfo->MaxNArr + iArr];
                    if(isAscii){
                        AARRFile << (float)factor * arr->a
                                 << (float)RadDeg * arr->Phase
                                 << arr->delay.real()
                                 << arr->delay.imag()
                                 << arr->SrcDeclAngle
                                 << arr->RcvrDeclAngle
                                 << arr->NTopBnc
                                 << arr->NBotBnc
                                 << '\n';
                    }else{
                        BARRFile.rec();
                        BARRFile.write((float)(factor * arr->a)); // LP: not a typo; cast in different order
                        BARRFile.write((float)(RadDeg * arr->Phase)); // compared to ascii version
                        BARRFile.write(arr->delay);
                        BARRFile.write(arr->SrcDeclAngle);
                        BARRFile.write(arr->RcvrDeclAngle);
                        BARRFile.write((float)arr->NTopBnc);
                        BARRFile.write((float)arr->NBotBnc);
                    }
                }
            }
        }
    }
}

}
