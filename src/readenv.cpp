#include "common.hpp"
#include "ldio.hpp"

#include <iostream>

/**
 * LP: No function description given.
 * 
 * bc: Boundary condition type
 */
void ReadTopOpt(char (TopOpt&)[6], char &bc, char (AttenUnit&)[2], 
    std::string FileRoot, LDIFile &ENVFile, std::ostream &PRTFile,
    SSPStructure &ssp, std::ifstream &SSPFile,
    AttenInfo &atten)
{
    TopOpt = "      "; // initialize to blanks
    ENVFile.Line(); ENVFile.Read(TopOpt, 6);
    PRTFile << "\n";
    
    ssp.Type     = TopOpt[0];
    bc           = TopOpt[1];
    AttenUnit[0] = TopOpt[2];
    AttenUnit[1] = TopOpt[3];
    ssp.AttenUnit = AttenUnit;
    
    // SSP approximation options
    
    switch(ssp.Type){
    case 'N':
        PRTFile << "    N2-linear approximation to SSP\n"; break;
    case 'C':
        PRTFile << "    C-linear approximation to SSP\n"; break;
    case 'P':
        PRTFile << "    PCHIP approximation to SSP\n"; break;
    case 'S':
        PRTFile << "    Spline approximation to SSP\n"; break;
    case 'Q':
        PRTFile << "    Quad approximation to SSP\n";
        SSPFile.open(WithExtension(FileRoot, ".ssp"));
        if(!SSPFile.good()){
            PRTFile << "SSPFile = " << WithExtension(FileRoot, ".ssp") << "\n";
            std::cout << "bellhopcuda - ReadEnvironment: Unable to open the SSP file\n";
        }
        break;
    /*case 'H':
        PRTFile << "    Hexahedral approximation to SSP\n";
        SSPFile.open(WithExtension(FileRoot, ".ssp"));
        if(!SSPFile.good()){
            PRTFile << "SSPFile = " << WithExtension(FileRoot, ".ssp") << "\n";
            std::cout << "bellhopcuda - ReadEnvironment: Unable to open the SSP file\n";
        }
        break;*/
    case 'A':
        PRTFile << "    Analytic SSP option\n"; break;
    default:
        std::cout << "ReadEnvironment: Unknown option for SSP approximation\n";
    }
    
    // Attenuation options
    
    switch(AttenUnit[0]){
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
    }
    
    // optional addition of volume attenuation using standard formulas
    
    switch(AttenUnit[1]){
    case 'T':
        PRTFile << "    THORP volume attenuation added\n"; break;
    case 'F':
        PRTFile << "    Francois-Garrison volume attenuation added\n";
        ENVFile.Line(); ENVFile.Read(atten.t); ENVFile.Read(atten.Salinity);
        ENVFile.Read(atten.pH); ENVFile.Read(atten.z_bar);
        PRTFile << std::setprecision(4);
        PRTFile << " T = " << std::setw(11) << atten.t 
            << "degrees   S = " << std::setw(11) << atten.Salinity
            << " psu   pH = " << std::setw(11) << atten.pH
            << " z_bar = " << std::setw(11) << " m\n";
        break;
    case 'B':
        PRTFile << "    Biological attenuation\n";
        ENVFile.Line(); ENVFile.Read(atten.NBioLayers);
        PRTFile << "      Number of Bio Layers = " << atten.NBioLayers << "\n";
        for(int32_t iBio = 0; iBio < NBioLayers; ++iBio){
            ENVFile.Line(); ENVFile.Read(atten.bio[iBio].z1); ENVFile.Read(atten.bio[iBio].z2);
            ENVFile.Read(atten.bio[iBio].f0); ENVFile.Read(atten.bio[iBio].q); ENVFile.Read(atten.bio[iBio].a0);
            PRTFile << "      Top    of layer = " << bio[iBio].z1 << " m\n";
            PRTFile << "      Bottom of layer = " << bio[iBio].z2 << " m\n";
            PRTFile << "      Resonance frequency = " << bio[iBio].f0 << " Hz\n";
            PRTFile << "      Q = " << bio[iBio].q << "\n";
            PRTFile << "      a0 = " << bio[iBio].a0 << "\n";
        }
    case ' ':
        break;
    default:
        std::cout << "ReadEnvironment: Unknown top option letter in fourth position\n";
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
    }
    
    switch(TopOpt[5]){
    case 'I':
        PRTFile << "    Development options enabled\n";
        break;
    case ' ':
        break;
    default:
        PRTFile << "ReadEnvironment: Unknown top option letter in sixth position\n";
    }
}

/**
 * Read the RunType variable and echo with explanatory information to the print file
 */
void ReadRunType(char (&RunType)[7], char (&PlotType)[10],
    LDIFile &ENVFile, std::ostream &PRTFile,
    Position &Pos)
{
    ENVFile.Line(); ENVFile.Read(RunType, 7);
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
       std::cout << "ReadEnvironment: Unknown RunType selected\n";
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
       if(Pos.NRz != Pos.NRr) std::cout << "ReadEnvironment: Irregular grid option selected with NRz not equal to Nr\n";
       PlotType = "irregular ";
       break;
    default:
       PRTFile << "Rectilinear receiver grid: Receivers at Rr[:] x Rz[:]\n";
       RunType[4] = 'R';
       PlotType = "rectilin  ";
    }

    switch(RunType[5]){
    case '2':
       PRTFile << "N x 2D calculation (neglects horizontal refraction)\n"; break;
    case '3':
       //PRTFile << "3D calculation\n";
       std::cout << "3D calculation not supported\n";
       RunType[5] = '2';
       break;
    default:
       RunType[5] = '2';
    }
}

/**
 * Handles top and bottom boundary conditions
 * 
 * freq: frequency [LP: :( ]
 */
void TopBot(const real &freq, const char (&AttenUnit)[2], HSInfo &hs,
    LDIFile &ENVFile, std::ostream &PRTFile,
    SSPVars &sspvars)
{
    real Mz, vr, alpha2_f; // values related to grain size
    
    // Echo to PRTFile user's choice of boundary condition
    
    switch(hs.bc){
    case 'V':
        PRTFile << "    VACUUM\n"; break;
    case 'R':
        PRTFile << "    Perfectly RIGID\n"; break;
    case 'A':
        PRTFile << "    ACOUSTO-ELASTIC half-space\n"; break;
    case 'G':
        PRTFile << "    Grain size to define half-space\n"; break;
    case 'F':
        PRTFile << "    FILE used for reflection loss\n"; break;
    case 'W':
        PRTFile << "    Writing an IRC file\n"; break;
    case 'P':
        PRTFile << "    reading PRECALCULATED IRC\n"; break;
    default:
       std::cout << "TopBot: Unknown boundary condition type\n";
    }
    
    // ****** Read in BC parameters depending on particular choice ******
    
    hs.cp = hs.cs = hs.rho = RC(0.0);
    
    if(hs.bc == 'A'){ // *** Half-space properties ***
        sspvars.zTemp = RC(0.0);
        ENVFile.List(); ENVFile.Read(sspvars.zTemp); ENVFile.Read(sspvars.alphaR);
        ENVFile.Read(sspvars.betaR); ENVFile.read(sspvars.rhoR);
        ENVFile.read(sspvars.alphaI); ENVFile.read(sspvars.beta);
        PRTFile << std::setprecision(2) << std::setw(10) << ztemp << " "
            << std::setw(10) << alphaR << " " << std::setw(10) << betaR << " "
            << std::setw(6) << rhoR << " " << std::setprecision(4) 
            << std::setw(10) << alphaI << " " << std::setw(10) << betaI << "\n";
        // dummy parameters for a layer with a general power law for attenuation
        // these are not in play because the AttenUnit for this is not allowed yet
        //freq0         = freq;
        sspvars.betaPowerLaw  = RC(1.0);
        sspvars.ft            = RC(1000.0);
        
        hs.cp  = crci(sspvars.zTemp, sspvars.alphaR, sspvars.alphaI, freq, freq,
            AttenUnit, sspvars.betaPowerLaw, sspvars.fT);
        hs.cs  = crci(sspvars.zTemp, sspvars.betaR,  sspvars.betaI,  freq, freq,
            AttenUnit, sspvars.betaPowerLaw, sspvars.fT);
        
        hs.rho = sspvars.rhoR;
    }else if(hs.bc == 'G'){ // *** Grain size (formulas from UW-APL HF Handbook)
        
        // These formulas are from the UW-APL Handbook
        // The code is taken from older Matlab and is unnecesarily verbose
        // vr   is the sound speed ratio
        // rhor is the density ratio
        ENVFile.List(); ENVFile.Read(sspvars.zTemp); ENVFile.Read(Mz);
        PRTFile << std::setprecision(2) << std::setw(10) << sspvars.zTemp << " "
            << std::setw(10) << Mz << "\n";
        
        if(Mz >= RC(-1.0) && Mz < RC(1.0)){
            vr           = RC(0.002709) * SQ(Mz) - RC(0.056452) * Mz + RC(1.2778);
            sspvars.rhor = RC(0.007797) * SQ(Mz) - RC(0.17057)  * Mz + RC(2.3139);
        }else if(Mz >= RC(1.0) && Mz < RC(5.3)){
            vr           = RC(-0.0014881) * Cube(Mz) + RC(0.0213937) * SQ(Mz) - RC(0.1382798) * Mz + RC(1.3425);
            sspvars.rhor = RC(-0.0165406) * Cube(Mz) + RC(0.2290201) * SQ(Mz) - RC(1.1069031) * Mz + RC(3.0455);
        }else{
            vr           = RC(-0.0024324) * Mz + RC(1.0019);
            sspvars.rhor = RC(-0.0012973) * Mz + RC(1.1565);
        }
        
        if(Mz >= RC(-1.0) && Mz < RC(0.0)){
            alpha2_f = RC(0.4556);
        }else if(Mz >= RC(0.0) && Mz < RC(2.6)){
            alpha2_f = RC(0.4556) + RC(0.0245) * Mz;
        }else if(Mz >= RC(2.6) && Mz < RC(4.5)){
            alpha2_f = RC(0.1978) + RC(0.1245) * Mz;
        }else if(Mz >= RC(4.5) && Mz < RC(6.0)){
            alpha2_f = RC(8.0399) - RC(2.5228) * Mz + RC(0.20098) * SQ(Mz);
        }else if(Mz >= RC(6.0) && Mz < RC(9.5)){
            alpha2_f = RC(0.9431) - RC(0.2041) * Mz + RC(0.0117) * SQ(Mz);
        }else{
            alpha2_f =  RC(0.0601);
        }
        
        // AttenUnit = 'L';  // loss parameter
        // !! following uses a reference sound speed of 1500 ???
        // !! should be sound speed in the water, just above the sediment
        // the term vr / 1000 converts vr to units of m per ms 
        sspvars.alphaR = vr * RC(1500.0);
        sspvars.alphaI = alpha2_f * (vr / RC(1000.0)) * RC(1500.0) * 
            STD::log(RC(10.0)) / (RC(40.0) * M_PI); // loss parameter Sect. IV., Eq. (4) of handbook
 
        hs.cp  = crci(sspvars.zTemp, sspvars.alphaR, sspvars.alphaI, freq, freq,
            "L ", sspvars.betaPowerLaw, sspvars.fT);
        hs.cs  = RC(0.0);
        hs.rho = sspvars.rhoR;
    }
}

void ReadEnvironment(std::string FileRoot, std::ostream &PRTFile,
    std::string &Title, real &freq, BdryType &Bdry, 
    SSPStructure &ssp, std::ifstream &SSPFile,
    AttenInfo &atten)
{
    const real c0 = RC(1500.0);
    int32_t NPts, NMedia;
    real ZMin, ZMax;
    vec2 x, gradc;
    cpx ccpx;
    real crr, crz, czz, rho, Sigma, Depth;
    char AttenUnit[2];
    char PlotType[10];
    
    PRTFile << "bellhopcuda\n\n";
    
    // Open the environmental file
    LDIFile ENVFile(WithExtension(FileRoot, ".env"));
    if(!ENVFile.Good()){
        PRTFile << "ENVFile = " << WithExtension(FileRoot, ".env") << "\n";
        std::cout << "bellhopcuda - ReadEnvironment: Unable to open the environmental file\n";
        return;
    }
    
    // Prepend model name to title
    ENVFile.List(); ENVFile.Read(Title);
    Title = "bellhopcuda- " + Title;
    
    PRTFile << Title << "\n";
    
    ENVFile.List(); ENVFile.Read(freq);
    PRTFile << std::setiosflags(std::scientific) << std::setprecision(4);
    PRTFile << " frequency = " << std::setw(11) << freq << " Hz\n");
    
    ENVFile.List(); ENVFile.Read(NMedia);
    PRTFile << "Dummy parameter NMedia = " << NMedia << "\n";
    if(NMedia != 1){
        std::cout << "ReadEnvironment: Only one medium or layer is allowed in BELLHOP; sediment layers must be handled using a reflection coefficient\n";
    }
    
    ReadTopOpt(Bdry.Top.hs.Opt, Bdr.Top.hs.bc, AttenUnit, FileRoot, 
        ENVFile, PRTFile, ssp, SSPFile, atten);
    
    // *** Top BC ***
    
    if(Bdry.Top.hs.bc == 'A') PRTFile << "   z (m)     alphaR (m/s)   betaR  rho (g/cm^3)  alphaI     betaI\n";
    
    TopBot(freq, AttenUnit, Bdry.Top.hs);
    
    // ****** Read in ocean SSP data ******
    
    ENVFile.List(); ENVFile.Read(NPts); ENVFile.Read(Sigma); ENVFile.Read(Bdry.Bot.hs.Depth);
    PRTFile << "\n  Depth = " << std::setw(10) << std::setprecision(2) << Bdr.Bot.hs.Depth << "  m )\n";
    
    if(Bdry.Top.hs.Opt[0] == 'A'){
        PRTFile << "Analytic SSP option\n";
        // following is hokey, should be set in Analytic routine
        ssp.Npts = 2;
        ssp.z.x = RC(0.0);
        ssp.z.y = Bdry.Bot.hs.Depth;
    }else{
        x = vec2(RC(0.0), Bdry.Bot.hs.Depth);
        InitializeSSP(x, freq);
    }
    
    Bdry.Top.hs.Depth = ssp.z[0]; // Depth of top boundary is taken from first SSP point
    // bottom depth should perhaps be set the same way?
    
    // *** Bottom BC ***
    Bdry.Bot.hs.Opt = "  "; // initialize to blanks
    ENVFile.List(); ENVFile.Read(Bdry.Bot.hs.Opt, 6); ENVFile.Read(Sigma);
    PRTFile << "\n RMS roughness = " << std::setw(10) << std::setprecision(3) << Sigma << "\n";
    
    switch(Bdry.Bot.hs.Opt[1]){
    case '~':
    case '*':
        PRTFile << "    Bathymetry file selected\n";
        break;
    case '-':
    case '_':
    case ' ':
        break;
    default:
        std::cout << "Unknown bottom option letter in second position\n";
    }
    
    Bdry.Bot.hs.bc = Bdr.Bot.hs.Opt[0];
    TopBot(freq, AttenUnit, Bdry.Bot.hs);
    
    // *** source and receiver locations ***
    
    ReadSxSy(false, ENVFile, PRTFile);
    
    ZMin = Bdry.Top.hs.Depth;
    ZMax = Bdry.Bot.hs.Depth;
    // not sure why I had this
    // ReadSzRz(ZMin + RC(100.0) * (std::nextafter(ZMin, ZMin+RC(1.0)) - ZMin),
    //     ZMax - RC(100.0) * (std::nextafter(ZMax, ZMax+RC(1.0)) - ZMax), 
    //     ENVFile, PRTFile);
    ReadSzRz(ZMin, ZMax, ENVFile, PRTFile);
    ReadRcvrRanges(ENVFile, PRTFile);
    ReadfreqVec(freq, Bdry.Top.hs.Opt[5], ENVFile, PRTFile);
    ReadRunType(Beam.RunType, PlotType, ENVFile, PRTFile);
    
    Depth = Zmax - Zmin; // water depth
    ReadRayElevationAngles(freq, Depth, Bdry.Top.hs.Opt, Beam.RunType, ENVFile, PRTFile);
    
    PRTFile << "\n__________________________________________________________________________\n\n";
    
    // Limits for tracing beams
    
    ENVFile.List(); ENVFile.Read(Beam.deltas); ENVFile.Read(Beam.Box.z); ENVFile.Read(Beam.Box.r);
    
    PRTFile << std::setprecision(4);
    PRTFile << "\n Step length,       deltas = " << std::setw(11) << Beam.deltas << " m\n\n";
    PRTFile << "Maximum ray depth, Box.z  = " << std::setw(11) << Beam.Box.z << " m\n";
    PRTFile << "Maximum ray range, Box.r  = " << std::setw(11) << Beam.Box.r << "km\n";
    
    Beam.Box.r *= RC(1000.0); // convert km to m
    
    // *** Beam characteristics ***
    
    Beam.Type[3] = Beam.RunType[6]; // selects beam shift option
    
    if(Beam.Type[3] == 'S'){
        PRTFile << "Beam shift in effect\n";
    }else{
        PRTFile << "No beam shift in effect\n";
    }
    
    if(Beam.RunType[0] != 'R'){ // no worry about the beam type if this is a ray trace run
        
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
        
        Beam.Type[0] = Beam.RunType[1];
        switch(Beam.Type[0]){
        case 'G':
        case 'g':
        case '^':
        case 'B':
        case 'b':
        case 'S':
            break;
        case 'R':
        case 'C':
            ENVFile.List(); ENVFile.Read(&Beam.Type[1], 2); ENVFile.Read(Beam.epsMultiplier); ENVFile.Read(Beam.rLoop);
            PRTFile << "\n\nType of beam = " << Beam.Type[0] << "\n";
            switch(Beam.Type[2]){
            case 'D':
                PRTFile << "Curvature doubling invoked\n"; break;
            case 'Z':
                PRTFile << "Curvature zeroing invoked\n"; break;
            case 'S':
                PRTFile << "Standard curvature condition\n"; break;
            default:
                std::cout << "ReadEnvironment: Unknown curvature condition\n";
            }
            break;
        default:
            std::cout << "ReadEnvironment: Unknown beam type (second letter of run type\n";
        }
    }
} 
