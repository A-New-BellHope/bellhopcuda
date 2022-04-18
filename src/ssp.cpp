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
#include "ssp.hpp"
#include "curves.hpp"
#include "boundary.hpp"

namespace bhc {

#define READ_SSP_ARGS real Depth, real freq, const real &fT, SSPStructure *ssp, \
    LDIFile &ENVFile, PrintFileEmu &PRTFile, const AttenInfo *atten, HSInfo &RecycledHS
#define CALL_READ_SSP_ARGS Depth, freqinfo->freq0, fT, ssp, ENVFile, PRTFile, atten, RecycledHS

/**
 * reads the SSP data from the environmental file and convert to Nepers/m
 */
void ReadSSP(READ_SSP_ARGS)
{
    PRTFile << "\nSound speed profile:\n";
    PRTFile << "   z (m)     alphaR (m/s)   betaR  rho (g/cm^3)  alphaI     betaI\n";
    
    ssp->NPts = 1;
    
    for(int32_t iz=0; iz<MaxSSP; ++iz){
        LIST_WARNLINE(ENVFile); ENVFile.Read(ssp->z[iz]); 
        ENVFile.Read(RecycledHS.alphaR); ENVFile.Read(RecycledHS.betaR);
        ENVFile.Read(RecycledHS.rho);
        ENVFile.Read(RecycledHS.alphaI); ENVFile.Read(RecycledHS.betaI);
        PRTFile << std::setprecision(2) << ssp->z[iz] << " " << RecycledHS.alphaR
            << " " << RecycledHS.betaR << " " << RecycledHS.rho << " " 
            << std::setprecision(4) << RecycledHS.alphaI << " " << RecycledHS.betaI << "\n";
        
        ssp->c[iz] = crci(ssp->z[iz], RecycledHS.alphaR, RecycledHS.alphaI, freq, freq,
            ssp->AttenUnit, betaPowerLaw, fT, atten, PRTFile);
        ssp->rho[iz] = RecycledHS.rho;
        
        // verify that the depths are monotone increasing
        if(iz > 0){
            if(ssp->z[iz] <= ssp->z[iz-1]){
                std::cout << "ReadSSP: The depths in the SSP must be monotone increasing (" << ssp->z[iz] << ")\n";
                std::abort();
            }
        }
        
        // compute gradient, cz
        if(iz > 0) ssp->cz[iz-1] = (ssp->c[iz] - ssp->c[iz-1]) /
                                   (ssp->z[iz] - ssp->z[iz-1]);
        
        // Did we read the last point?
        if(std::abs(ssp->z[iz] - Depth) < FL(100.0) * FLT_EPSILON){ // LP: FLT_EPSILON is not a typo
            // LP: Gradient at iz is uninitialized.
            ssp->cz[iz] = cpx(FL(5.5555555e30), FL(-3.3333333e29)); // LP: debugging
            
            ssp->Nz = ssp->NPts;
            if(ssp->NPts == 1){
                std::cout << "ReadSSP: The SSP must have at least 2 points\n";
                std::abort();
            }
            
            return;
        }
        
        ++ssp->NPts;
    }
    
    // Fall through means too many points in the profile
    std::cout << "ReadSSP: Number of SSP points exceeds limit\n";
    std::abort();
}

void Initn2Linear(SSP_INIT_ARGS)
{
    IGNORE_UNUSED(FileRoot);
    
    ReadSSP(CALL_READ_SSP_ARGS);
    
    for(int32_t i=0; i<ssp->NPts; ++i) ssp->n2[i] = FL(1.0) / SQ(ssp->c[i]);
    
    // compute gradient, n2z
    for(int32_t iz=1; iz<ssp->NPts; ++iz){
        ssp->n2z[iz-1] = (ssp->n2[iz] - ssp->n2[iz-1]) /
                         (ssp->z [iz] - ssp->z [iz-1]);
    }
}

void InitcLinear(SSP_INIT_ARGS)
{
    IGNORE_UNUSED(FileRoot);
    
    ReadSSP(CALL_READ_SSP_ARGS);
}

void InitcPCHIP(SSP_INIT_ARGS)
{
    IGNORE_UNUSED(FileRoot);
    
    ReadSSP(CALL_READ_SSP_ARGS);
    
    //                                                               2      3
    // compute coefficients of std cubic polynomial: c0 + c1*x + c2*x + c3*x
    //
    pchip(ssp->z, ssp->c, ssp->NPts, 
        ssp->cCoef[0], ssp->cCoef[1], ssp->cCoef[2], ssp->cCoef[3],
        ssp->CSWork[0], ssp->CSWork[1], ssp->CSWork[2], ssp->CSWork[3]);
}

void InitcCubic(SSP_INIT_ARGS)
{
    IGNORE_UNUSED(FileRoot);
    
    ReadSSP(CALL_READ_SSP_ARGS);
    
    for(int32_t i=0; i<ssp->NPts; ++i) ssp->cSpline[0][i] = ssp->c[i];
    
    // Compute spline coefs
    int32_t iBCBeg = 0;
    int32_t iBCEnd = 0;
    cSpline(ssp->z, ssp->cSpline[0], ssp->cSpline[1], ssp->cSpline[2], ssp->cSpline[3],
        ssp->NPts, iBCBeg, iBCEnd, ssp->NPts);
}

void InitQuad(SSP_INIT_ARGS)
{
    ReadSSP(CALL_READ_SSP_ARGS);
    
    // Read the 2D SSP matrix
    PRTFile << "__________________________________________________________________________\n\n";
    PRTFile << "Using range-dependent sound speed\n";
    
    LDIFile SSPFile(FileRoot + ".ssp");
    LIST(SSPFile); SSPFile.Read(ssp->Nr);
    PRTFile << "Number of SSP ranges = " << ssp->Nr << "\n";
    
    if(ssp->Nr < 2){
        PRTFile << "READIN: Quad: You must have a least two profiles in your 2D SSP field\n";
        std::abort();
    }
    
    checkallocate(ssp->cMat ,  ssp->NPts    * ssp->Nr);
    checkallocate(ssp->czMat, (ssp->NPts-1) * ssp->Nr);
    checkallocate(ssp->Seg.r, ssp->Nr);
    
    LIST(SSPFile); SSPFile.Read(ssp->Seg.r, ssp->Nr);
    PRTFile << "\nProfile ranges (km):\n" << std::setprecision(2);
    for(int32_t i=0; i<ssp->Nr; ++i) PRTFile << ssp->Seg.r[i] << " ";
    PRTFile << "\n";
    
    for(int32_t i=0; i<ssp->Nr; ++i) ssp->Seg.r[i] *= FL(1000.0); // convert km to m
    
    PRTFile << "\nSound speed matrix:\n";
    PRTFile << " Depth (m )     Soundspeed (m/s)\n";
    for(int32_t iz2=0; iz2<ssp->NPts; ++iz2){
        LIST(SSPFile); SSPFile.Read(&ssp->cMat[iz2*ssp->Nr], ssp->Nr);
        // PRTFile << "iSeg.z depth = " << std::setprecision(2) << ssp->z[iz2] << " m\n";
        PRTFile << std::setprecision(2) << ssp->z[iz2] << " ";
        for(int32_t i=0; i<ssp->Nr; ++i) PRTFile << ssp->cMat[iz2*ssp->Nr+i] << " ";
        PRTFile << "\n";
    }
    
    // calculate cz
    for(int32_t iSegt=0; iSegt<ssp->Nr; ++iSegt){
        for(int32_t iz2=1; iz2<ssp->NPts; ++iz2){
            real delta_z = ssp->z[iz2] - ssp->z[iz2-1];
            ssp->czMat[(iz2-1)*ssp->Nr + iSegt] = 
                (ssp->cMat[iz2*ssp->Nr+iSegt] - ssp->cMat[(iz2-1)*ssp->Nr+iSegt]) / delta_z;
        }
    }
    
    ssp->Nz = ssp->NPts;
}

void InitHexahedral(SSP_INIT_ARGS)
{
    // Read dummy SSP information from the environmental file
    // This is over-ridden by the info in the SSP file
    // However, cz info is used in beam selection
    
    ReadSSP(CALL_READ_SSP_ARGS);
    
    // Read the 3D SSP matrix
    
    PRTFile << "\nReading sound speed profile from file\n";
    LDIFile SSPFile(FileRoot + ".ssp");
    
    // x coordinates
    LIST(SSPFile); SSPFile.Read(ssp->Nx);
    PRTFile << "\nNumber of points in x = " << ssp->Nx << "\n";
    checkallocate(ssp->Seg.x, ssp->Nx);
    LIST(SSPFile); SSPFile.read(ssp->Seg.x, ssp->Nx);
    
    // y coordinates
    LIST(SSPFile); SSPFile.Read(ssp->Ny);
    PRTFile << "\nNumber of points in y = " << ssp->Ny << "\n";
    checkallocate(ssp->Seg.y, ssp->Ny);
    LIST(SSPFile); SSPFile.read(ssp->Seg.y, ssp->Ny);
    
    // z coordinates
    LIST(SSPFile); SSPFile.Read(ssp->Nz);
    PRTFile << "\nNumber of points in z = " << ssp->Nz << "\n";
    checkallocate(ssp->Seg.z, ssp->Nz);
    LIST(SSPFile); SSPFile.read(ssp->Seg.z, ssp->Nz);
    
    if(ssp->Nx < 2 || ssp->Ny < 2 || ssp->Nz < 2){
        std::cout << "You must have at least two points in x, y, z directions in your 3D SSP field\n";
        std::abort();
    }
    
    checkallocate(ssp->cMat,  ssp->Nx * ssp->Ny * ssp->Nz);
    checkallocate(ssp->czMat, ssp->Nx * ssp->Ny * (ssp->Nz - 1));
    
    PRTFile << "\n";
    for(int32_t iz2=0; iz2<ssp->Nz; ++iz2){
        for(int32_t iy2=0; iy2<ssp->Ny; ++iy2){
            LIST(SSPFile);
            for(int32_t ix2=0; ix2<ssp->Nx; ++ix2){
                SSPFile.read(&ssp->cMat[((ix2)*ssp->Ny+iy2)*ssp->Nz+iz2]);
            }
        }
    }
    
    // convert km to m
    for(int32_t ix1=0; ix1<ssp->Nx; ++ix1) ssp->Seg.x *= FL(1000.0);
    for(int32_t iy1=0; iy1<ssp->Ny; ++iy1) ssp->Seg.y *= FL(1000.0);
    
    // calculate cz
    for(int32_t iSegxt=0; iSegxt<ssp->Nx; ++iSegxt){
        for(int32_t iSegyt=0; iSegyt<ssp->Ny; ++iSegyt){
            for(int32_t iz2=1; iz2<ssp->Nz; ++iz2){
                ssp->czMat[((iSegxt)*ssp->Ny+iSegyt)*(ssp->Nz-1)+iz2-1] =
                    (ssp->cMat[((iSegxt)*ssp->Ny+iSegyt)*ssp->Nz+iz2] - ssp->cMat[((iSegxt)*ssp->Ny+iSegyt)*ssp->Nz+iz2-1]) /
                    (ssp->Seg.z[                                 iz2] - ssp->Seg.z[                                 iz2-1]);
            }
        }
    }
    
    // over-ride the SSP%z values read in from the environmental file with these new values
    ssp->NPts = ssp->Nz;
    for(int32_t iz=0; iz<ssp->Nz; ++iz) ssp->z[iz] = ssp->Seg.z[iz];
}

void InitializeSSP(SSP_INIT_ARGS)
{
    switch(ssp->Type){
    case 'N': // N2-linear profile option
        Initn2Linear  (SSP_CALL_INIT_ARGS); break;
    case 'C': // C-linear profile option
        InitcLinear   (SSP_CALL_INIT_ARGS); break;
    case 'P': // monotone PCHIP ACS profile option
        InitcPCHIP    (SSP_CALL_INIT_ARGS); break;
    case 'S': // Cubic spline profile option
        InitcCubic    (SSP_CALL_INIT_ARGS); break;
    case 'Q':
        InitQuad      (SSP_CALL_INIT_ARGS); break;
    case 'H':
        InitHexahedral(SSP_CALL_INIT_ARGS); break;
    case 'A': // Analytic profile option
        break; //LP: No init for analytic 2D or 3D.
    default:
        printf("InitializeSSP: Invalid profile option %c\n", ssp->Type);
        std::abort();
    }
}

}
