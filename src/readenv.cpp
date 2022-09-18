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
#include "readenv.hpp"
#include "boundary.hpp"
#include "ssp.hpp"
#include "sourcereceiver.hpp"
#include "angles.hpp"
#include "beams.hpp"

namespace bhc {

/**
 * LP: No function description given.
 * 
 * bc: Boundary condition type
 */
void ReadTopOpt(char (&TopOpt)[6], char &bc,
    std::string FileRoot, LDIFile &ENVFile, PrintFileEmu &PRTFile,
    SSPStructure *ssp, AttenInfo *atten)
{
    memcpy(TopOpt, "      ", 6); // initialize to blanks
    LIST(ENVFile); ENVFile.Read(TopOpt, 6);
    PRTFile << "\n";
    
    ssp->Type         = TopOpt[0];
    bc                = TopOpt[1];
    ssp->AttenUnit[0] = TopOpt[2];
    ssp->AttenUnit[1] = TopOpt[3];
    
    // SSP approximation options
    
    switch(ssp->Type){
    case 'N':
        PRTFile << "    N2-linear approximation to SSP\n"; break;
    case 'C':
        PRTFile << "    C-linear approximation to SSP\n"; break;
    case 'P':
        PRTFile << "    PCHIP approximation to SSP\n"; break;
    case 'S':
        PRTFile << "    Spline approximation to SSP\n"; break;
    case 'Q':{
        PRTFile << "    Quad approximation to SSP\n";
        //LP: This just checks for existence, moved actual open for reading
        //to InitQuad.
        std::ifstream SSPFile;
        SSPFile.open(FileRoot + ".ssp");
        if(!SSPFile.good()){
            PRTFile << "SSPFile = " << FileRoot << ".ssp\n";
            std::cout << BHC_PROGRAMNAME " - ReadEnvironment: Unable to open the SSP file\n";
            std::abort();
        }
        } break;
    case 'H':{
        //LP: This just checks for existence, moved actual open for reading
        //to InitHexahedral.
        PRTFile << "    Hexahedral approximation to SSP\n";
        std::ifstream SSPFile;
        SSPFile.open(FileRoot + ".ssp");
        if(!SSPFile.good()){
            PRTFile << "SSPFile = " << FileRoot << ".ssp\n";
            std::cout << BHC_PROGRAMNAME " - ReadEnvironment: Unable to open the SSP file\n";
        }
        } break;
    case 'A':
        PRTFile << "    Analytic SSP option\n"; break;
    default:
        std::cout << "ReadEnvironment: Unknown option for SSP approximation\n";
        std::abort();
    }
    
    // Attenuation options
    
    switch(ssp->AttenUnit[0]){
    case 'N':
        PRTFile << "    Attenuation units: nepers/m\n"; break;
    case 'F':
        PRTFile << "    Attenuation units: dB/mkHz\n"; break;
    case 'M':
        PRTFile << "    Attenuation units: dB/m\n"; break;
    case 'W':
        PRTFile << "    Attenuation units: dB/wavelength\n"; break;
    case 'Q':
        PRTFile << "    Attenuation units: Q\n"; break;
    case 'L':
        PRTFile << "    Attenuation units: Loss parameter\n"; break;
    default:
        std::cout << "ReadEnvironment: Unknown attenuation units\n";
        std::abort();
    }
    
    // optional addition of volume attenuation using standard formulas
    
    switch(ssp->AttenUnit[1]){
    case 'T':
        PRTFile << "    THORP volume attenuation added\n"; break;
    case 'F':
        PRTFile << "    Francois-Garrison volume attenuation added\n";
        LIST(ENVFile); ENVFile.Read(atten->t); ENVFile.Read(atten->Salinity);
        ENVFile.Read(atten->pH); ENVFile.Read(atten->z_bar);
        PRTFile << std::setprecision(4);
        PRTFile << " T = " << std::setw(11) << atten->t 
            << "degrees   S = " << std::setw(11) << atten->Salinity
            << " psu   pH = " << std::setw(11) << atten->pH
            << " z_bar = " << std::setw(11) << " m\n";
        break;
    case 'B':
        PRTFile << "    Biological attenuation\n";
        LIST(ENVFile); ENVFile.Read(atten->NBioLayers);
        PRTFile << "      Number of Bio Layers = " << atten->NBioLayers << "\n";
        for(int32_t iBio = 0; iBio < atten->NBioLayers; ++iBio){
            LIST(ENVFile); ENVFile.Read(atten->bio[iBio].z1); ENVFile.Read(atten->bio[iBio].z2);
            ENVFile.Read(atten->bio[iBio].f0); ENVFile.Read(atten->bio[iBio].q); ENVFile.Read(atten->bio[iBio].a0);
            PRTFile << "      Top    of layer = " << atten->bio[iBio].z1 << " m\n";
            PRTFile << "      Bottom of layer = " << atten->bio[iBio].z2 << " m\n";
            PRTFile << "      Resonance frequency = " << atten->bio[iBio].f0 << " Hz\n";
            PRTFile << "      Q = " << atten->bio[iBio].q << "\n";
            PRTFile << "      a0 = " << atten->bio[iBio].a0 << "\n";
        }
    case ' ':
        break;
    default:
        std::cout << "ReadEnvironment: Unknown top option letter in fourth position\n";
        std::abort();
    }
    
    switch(TopOpt[4]){
    case '~':
    case '*':
        PRTFile << "    Altimetry file selected\n";
        break;
    case '-':
    case '_':
    case ' ':
        break;
    default:
        PRTFile << "ReadEnvironment: Unknown top option letter in fifth position\n";
        std::abort();
    }
    
    switch(TopOpt[5]){
    case 'I':
        PRTFile << "    Development options enabled\n";
        break;
    case ' ':
        break;
    default:
        PRTFile << "ReadEnvironment: Unknown top option letter in sixth position\n";
        std::abort();
    }
}

/**
 * Read the RunType variable and echo with explanatory information to the print file
 */
void ReadRunType(char (&RunType)[7], char (&PlotType)[10],
    LDIFile &ENVFile, PrintFileEmu &PRTFile,
    Position *Pos, bool r3d)
{
    LIST(ENVFile); ENVFile.Read(RunType, 7);
    PRTFile << "\n";
    
    switch(RunType[0]){
    case 'R':
       PRTFile << "Ray trace run\n"; break;
    case 'E':
       PRTFile << "Eigenray trace run\n"; break;
    case 'I':
       PRTFile << "Incoherent TL calculation\n"; break;
    case 'S':
       PRTFile << "Semi-coherent TL calculation\n"; break;
    case 'C':
       PRTFile << "Coherent TL calculation\n"; break;
    case 'A':
       PRTFile << "Arrivals calculation, ASCII  file output\n"; break;
    case 'a':
       PRTFile << "Arrivals calculation, binary file output\n"; break;
    default:
       std::cout << "ReadEnvironment: Unknown RunType '" << RunType[0] << "' selected\n";
       std::abort();
    }

    switch(RunType[1]){
    case 'C':
       PRTFile << "Cartesian beams\n"; break;
    case 'R':
       PRTFile << "Ray centered beams\n"; break;
    case 'S':
       PRTFile << "Simple gaussian beams\n"; break;
    case 'b':
       PRTFile << "Geometric gaussian beams in ray-centered coordinates\n"; break;
    case 'B':
       PRTFile << "Geometric gaussian beams in Cartesian coordinates\n"; break;
    case 'g':
       PRTFile << "Geometric hat beams in ray-centered coordinates\n"; break;
    default:
       RunType[1] = 'G';
       PRTFile << "Geometric hat beams in Cartesian coordinates\n";
    }

    switch(RunType[3]){
    case 'X':
       PRTFile << "Line source (Cartesian coordinates)\n"; break;
    default:
       RunType[3] = 'R';
       PRTFile << "Point source (cylindrical coordinates)\n";
    }

    switch(RunType[4]){
    case 'I':
       PRTFile << "Irregular grid: Receivers at Rr[:] x Rz[:]\n";
       if(Pos->NRz != Pos->NRr) std::cout << "ReadEnvironment: Irregular grid option selected with NRz not equal to Nr\n";
       memcpy(PlotType, "irregular ", 10);
       break;
    default:
       PRTFile << "Rectilinear receiver grid: Receivers at Rr[:] x Rz[:]\n";
       RunType[4] = 'R';
       memcpy(PlotType, "rectilin  ", 10);
    }

    switch(RunType[5]){
    case '2':
       PRTFile << "N x 2D calculation (neglects horizontal refraction)\n";
       if(r3d){
           std::cout << "This is a 2D or Nx2D environment file, but you are running " BHC_PROGRAMNAME " in 3D mode\n";
           std::abort();
       }
       break;
    case '3':
       PRTFile << "3D calculation\n";
       if(!r3d){
           std::cout << "This is a 3D environment file, but you are running " BHC_PROGRAMNAME " in 2D or Nx2D mode\n";
           std::abort();
       }
       break;
    default:
       RunType[5] = '2';
    }
}

void ReadEnvironment(const std::string &FileRoot, PrintFileEmu &PRTFile,
    char (&Title)[80], real &fT, BdryType *Bdry, SSPStructure *ssp, AttenInfo *atten, 
    Position *Pos, AnglesStructure *Angles, FreqInfo *freqinfo, BeamStructure *Beam,
    HSInfo &RecycledHS, bool o3d, bool r3d)
{
    //const real c0 = FL(1500.0); //LP: unused
    int32_t NPts, NMedia;
    real ZMin, ZMax;
    vec2 x;
    cpx ccpx;
    real Sigma, Depth;
    char PlotType[10];
    
    PRTFile << BHC_PROGRAMNAME << (r3d ? "3D" : o3d ? "Nx2D" : "") << "\n\n";
    
    // Open the environmental file
    LDIFile ENVFile(FileRoot + ".env");
    if(!ENVFile.Good()){
        PRTFile << "ENVFile = " << FileRoot << ".env\n";
        std::cout << BHC_PROGRAMNAME " - ReadEnvironment: Unable to open the environmental file\n";
        std::abort();
    }
    
    // Prepend model name to title
    std::string TempTitle;
    LIST(ENVFile); ENVFile.Read(TempTitle);
    TempTitle = BHC_PROGRAMNAME "- " + TempTitle;
    PRTFile << TempTitle << "\n";
    size_t l = bhc::min(sizeof(Title) - 1, TempTitle.size());
    memcpy(Title, TempTitle.c_str(), l);
    Title[l] = 0;
    
    LIST(ENVFile); ENVFile.Read(freqinfo->freq0);
    PRTFile << std::setiosflags(std::ios::scientific) << std::setprecision(4);
    PRTFile << " frequency = " << std::setw(11) << freqinfo->freq0 << " Hz\n";
    
    LIST(ENVFile); ENVFile.Read(NMedia);
    PRTFile << "Dummy parameter NMedia = " << NMedia << "\n";
    if(NMedia != 1){
        std::cout << "ReadEnvironment: Only one medium or layer is allowed in BELLHOP; sediment layers must be handled using a reflection coefficient\n";
    }
    
    ReadTopOpt(Bdry->Top.hs.Opt, Bdry->Top.hs.bc, FileRoot, 
        ENVFile, PRTFile, ssp, atten);
    
    // *** Top BC ***
    
    if(Bdry->Top.hs.bc == 'A') PRTFile << "   z (m)     alphaR (m/s)   betaR  rho (g/cm^3)  alphaI     betaI\n";
    
    TopBot(freqinfo->freq0, ssp->AttenUnit, fT, Bdry->Top.hs, ENVFile, PRTFile, atten, RecycledHS);
    
    // ****** Read in ocean SSP data ******
    
    LIST(ENVFile); ENVFile.Read(NPts); ENVFile.Read(Sigma); ENVFile.Read(Bdry->Bot.hs.Depth);
    PRTFile << "\n  Depth = " << std::setw(10) << std::setprecision(2) << Bdry->Bot.hs.Depth << "  m\n";
    
    if(Bdry->Top.hs.Opt[0] == 'A'){
        PRTFile << "Analytic SSP option\n";
        // following is hokey, should be set in Analytic routine
        ssp->NPts = 2;
        ssp->z[0] = FL(0.0);
        ssp->z[1] = Bdry->Bot.hs.Depth;
    }else{
        real Depth = Bdry->Bot.hs.Depth;
        InitializeSSP(SSP_CALL_INIT_ARGS);
    }
    
    Bdry->Top.hs.Depth = ssp->z[0]; // Depth of top boundary is taken from first SSP point
    // bottom depth should perhaps be set the same way?
    
    // *** Bottom BC ***
    memcpy(Bdry->Bot.hs.Opt, "      ", 6); // initialize to blanks
    LIST(ENVFile); ENVFile.Read(Bdry->Bot.hs.Opt, 6); ENVFile.Read(Sigma);
    PRTFile << "\n RMS roughness = " << std::setw(10) << std::setprecision(3) << Sigma << "\n";
    
    switch(Bdry->Bot.hs.Opt[1]){
    case '~':
    case '*':
        PRTFile << "    Bathymetry file selected\n";
        break;
    case '-':
    case '_':
    case ' ':
        break;
    default:
        std::cout << "Unknown bottom option letter in second position: Bdr->Bot.hs.Opt == '" 
            << Bdry->Bot.hs.Opt << "'\n";
        std::abort();
    }
    
    Bdry->Bot.hs.bc = Bdry->Bot.hs.Opt[0];
    TopBot(freqinfo->freq0, ssp->AttenUnit, fT, Bdry->Bot.hs, ENVFile, PRTFile, atten, RecycledHS);
    
    // *** source and receiver locations ***
    
    ReadSxSy(o3d, ENVFile, PRTFile, Pos);
    
    ZMin = Bdry->Top.hs.Depth;
    ZMax = Bdry->Bot.hs.Depth;
    // not sure why I had this
    // ReadSzRz(ZMin + FL(100.0) * spacing(ZMin),
    //     ZMax - FL(100.0) * spacing(ZMax), 
    //     ENVFile, PRTFile);
    ReadSzRz(ZMin, ZMax, ENVFile, PRTFile, Pos);
    ReadRcvrRanges(ENVFile, PRTFile, Pos);
    if(o3d) ReadRcvrBearings(ENVFile, PRTFile, Pos);
    ReadfreqVec(Bdry->Top.hs.Opt[5], ENVFile, PRTFile, freqinfo);
    ReadRunType(Beam->RunType, PlotType, ENVFile, PRTFile, Pos, r3d);
    
    Depth = ZMax - ZMin; // water depth
    ReadRayAngles<false>(freqinfo->freq0, Depth, Bdry->Top.hs.Opt, 
        ENVFile, PRTFile, Angles->alpha, Pos, Beam);
    if(o3d) ReadRayAngles<true>(freqinfo->freq0, Depth, Bdry->Top.hs.Opt, 
        ENVFile, PRTFile, Angles->beta, Pos, Beam);
    
    PRTFile << "\n__________________________________________________________________________\n\n";
        
    //LP: Moved to separate function for clarity and modularity.
    ReadBeamInfo(ENVFile, PRTFile, Beam, Bdry, o3d);
} 

}
