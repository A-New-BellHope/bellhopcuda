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
#include "tlmode.hpp"

namespace bhc {

/**
 * Write header to disk file
 * LP: of a SHDFile ("binary `shade' file (SHDFIL) that contains calculated
 * pressure fields")
 *
 * FileName: Name of the file (could be a shade file or a Green's function file)
 * Title: Arbitrary title
 * freq0: Nominal frequency [LP: now in freqinfo]
 * Atten: stabilizing attenuation (for wavenumber integration only)
 * PlotType: If "TL", writes only first and last Sx and Sy [LP: never set to
 * "TL" in BELLHOP]
 */
void WriteHeader(DirectOFile &SHDFile, const std::string &FileName, 
    const char (&Title)[80], float Atten, const std::string &PlotType, 
    const Position *Pos, const FreqInfo *freqinfo)
{
    bool isTL = (PlotType[0] == 'T' && PlotType[1] == 'L');
    
    int32_t LRecl = 84; //4 for LRecl, 80 for Title
    LRecl = bhc::max(LRecl, 2 * freqinfo->Nfreq * (int32_t)sizeof(freqinfo->freqVec[0]));
    LRecl = bhc::max(LRecl, Pos->Ntheta * (int32_t)sizeof(Pos->theta[0]));
    if(!isTL){
        LRecl = bhc::max(LRecl, Pos->NSx * (int32_t)sizeof(Pos->Sx[0]));
        LRecl = bhc::max(LRecl, Pos->NSy * (int32_t)sizeof(Pos->Sy[0]));
    }
    LRecl = bhc::max(LRecl, Pos->NSz * (int32_t)sizeof(Pos->Sz[0]));
    LRecl = bhc::max(LRecl, Pos->NRz * (int32_t)sizeof(Pos->Rz[0]));
    LRecl = bhc::max(LRecl, Pos->NRr * (int32_t)sizeof(cpxf));
    
    SHDFile.open(FileName, LRecl);
    if(!SHDFile.good()){
        std::cout << "Could not open SHDFile: " << FileName << "\n";
        std::abort();
    }
    LRecl /= 4;
    SHDFile.rec(0); DOFWRITE(SHDFile, &LRecl, 4); DOFWRITE(SHDFile, std::string(Title), 80);
    SHDFile.rec(1); DOFWRITE(SHDFile, PlotType, 10);
    SHDFile.rec(2);
        DOFWRITEV(SHDFile, freqinfo->Nfreq);
        DOFWRITEV(SHDFile, Pos->Ntheta);
        DOFWRITEV(SHDFile, Pos->NSx);
        DOFWRITEV(SHDFile, Pos->NSy);
        DOFWRITEV(SHDFile, Pos->NSz);
        DOFWRITEV(SHDFile, Pos->NRz);
        DOFWRITEV(SHDFile, Pos->NRr);
        DOFWRITEV(SHDFile, (float)freqinfo->freq0);
        DOFWRITEV(SHDFile, Atten);
    SHDFile.rec(3); DOFWRITE(SHDFile, freqinfo->freqVec, freqinfo->Nfreq * sizeof(freqinfo->freqVec[0]));
    SHDFile.rec(4); DOFWRITE(SHDFile, Pos->theta, Pos->Ntheta * sizeof(Pos->theta[0]));
    
    if(!isTL){
        SHDFile.rec(5); DOFWRITE(SHDFile, Pos->Sx, Pos->NSx * sizeof(Pos->Sx[0]));
        SHDFile.rec(6); DOFWRITE(SHDFile, Pos->Sy, Pos->NSy * sizeof(Pos->Sy[0]));
    }else{
        SHDFile.rec(5); DOFWRITEV(SHDFile, Pos->Sx[0]); DOFWRITEV(SHDFile, Pos->Sx[Pos->NSx-1]);
        SHDFile.rec(6); DOFWRITEV(SHDFile, Pos->Sy[0]); DOFWRITEV(SHDFile, Pos->Sy[Pos->NSy-1]);
    }
    SHDFile.rec(7); DOFWRITE(SHDFile, Pos->Sz, Pos->NSz * sizeof(Pos->Sz[0]));
    
    SHDFile.rec(8); DOFWRITE(SHDFile, Pos->Rz, Pos->NRz * sizeof(Pos->Rz[0]));
    SHDFile.rec(9); DOFWRITE(SHDFile, Pos->Rr, Pos->NRr * sizeof(Pos->Rr[0]));
}


/**
 * LP: Write TL results
 */
void FinalizeTLMode(std::string FileRoot, const bhcParams &params, bhcOutputs &outputs)
{
    real atten = FL(0.0);
    std::string PlotType;
    DirectOFile SHDFile;
    
    // following to set PlotType has already been done in READIN if that was used for input
    switch(params.Beam->RunType[4]){
    case 'R':
        PlotType = "rectilin  "; break;
    case 'I':
        PlotType = "irregular "; break;
    default:
        PlotType = "rectilin  ";
    }
    WriteHeader(SHDFile, FileRoot + ".shd", params.Title, atten, PlotType,
        params.Pos, params.freqinfo);
    
    for(int32_t isrc=0; isrc<params.Pos->NSz; ++isrc){
        SSPSegState iSeg; iSeg.r = 0; iSeg.z = 0;
        SSPOutputs<false> o;
        Origin<false, false> org;
        EvaluateSSP<false, false>(vec2(RL(0.0), params.Pos->Sz[isrc]), vec2(RL(1.0), RL(0.0)),
            o, org, params.ssp, iSeg);
        ScalePressure(params.Angles->alpha.d, o.ccpx.real(), params.Pos->Rr, 
            &outputs.uAllSources[isrc * params.Pos->NRz_per_range * params.Pos->NRr], 
            params.Pos->NRz_per_range, params.Pos->NRr, params.Beam->RunType,
            params.freqinfo->freq0);
        int32_t IRec = 10 + params.Pos->NRz_per_range * isrc;
        for(int32_t Irz1 = 0; Irz1 < params.Pos->NRz_per_range; ++Irz1){
            SHDFile.rec(IRec);
            for(int32_t r=0; r < params.Pos->NRr; ++r){
                DOFWRITEV(SHDFile, outputs.uAllSources[
                    (isrc * params.Pos->NRz_per_range + Irz1) * params.Pos->NRr + r]);
            }
            ++IRec;
        }
    }
    
    deallocate(outputs.uAllSources);
}

}
