#include "readenv.hpp"
#include "ldio.hpp"

/**
 * LP: No function description given.
 * 
 * bc: Boundary condition type
 */
void ReadTopOpt(char (&TopOpt)[6], char &bc,
    std::string FileRoot, LDIFile &ENVFile, std::ofstream &PRTFile,
    SSPStructure *ssp, AttenInfo *atten)
{
    memcpy(TopOpt, "      ", 6); // initialize to blanks
    ENVFile.List(); ENVFile.Read(TopOpt, 6);
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
            std::cout << "bellhopcuda - ReadEnvironment: Unable to open the SSP file\n";
        }
        } break;
    /*case 'H':{
        PRTFile << "    Hexahedral approximation to SSP\n";
        std::ifstream SSPFile;
        SSPFile.open(FileRoot + ".ssp");
        if(!SSPFile.good()){
            PRTFile << "SSPFile = " << FileRoot << ".ssp\n";
            std::cout << "bellhopcuda - ReadEnvironment: Unable to open the SSP file\n";
        }
        } break;*/
    case 'A':
        PRTFile << "    Analytic SSP option\n"; break;
    default:
        std::cout << "ReadEnvironment: Unknown option for SSP approximation\n";
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
    }
    
    // optional addition of volume attenuation using standard formulas
    
    switch(ssp->AttenUnit[1]){
    case 'T':
        PRTFile << "    THORP volume attenuation added\n"; break;
    case 'F':
        PRTFile << "    Francois-Garrison volume attenuation added\n";
        ENVFile.List(); ENVFile.Read(atten->t); ENVFile.Read(atten->Salinity);
        ENVFile.Read(atten->pH); ENVFile.Read(atten->z_bar);
        PRTFile << std::setprecision(4);
        PRTFile << " T = " << std::setw(11) << atten->t 
            << "degrees   S = " << std::setw(11) << atten->Salinity
            << " psu   pH = " << std::setw(11) << atten->pH
            << " z_bar = " << std::setw(11) << " m\n";
        break;
    case 'B':
        PRTFile << "    Biological attenuation\n";
        ENVFile.List(); ENVFile.Read(atten->NBioLayers);
        PRTFile << "      Number of Bio Layers = " << atten->NBioLayers << "\n";
        for(int32_t iBio = 0; iBio < atten->NBioLayers; ++iBio){
            ENVFile.List(); ENVFile.Read(atten->bio[iBio].z1); ENVFile.Read(atten->bio[iBio].z2);
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
    LDIFile &ENVFile, std::ofstream &PRTFile,
    Position *Pos)
{
    ENVFile.List(); ENVFile.Read(RunType, 7);
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

void ReadEnvironment(const std::string &FileRoot, std::ofstream &PRTFile,
    std::string &Title, real &fT, BdryType *Bdry, SSPStructure *ssp, AttenInfo *atten, 
    Position *Pos, AnglesStructure *Angles, FreqInfo *freqinfo, BeamStructure *Beam)
{
    //const real c0 = RC(1500.0); //LP: unused
    int32_t NPts, NMedia;
    real ZMin, ZMax;
    vec2 x, gradc;
    cpx ccpx;
    real Sigma, Depth;
    char PlotType[10];
    
    PRTFile << "bellhopcuda\n\n";
    
    // Open the environmental file
    LDIFile ENVFile(FileRoot + ".env");
    if(!ENVFile.Good()){
        PRTFile << "ENVFile = " << FileRoot << ".env\n";
        std::cout << "bellhopcuda - ReadEnvironment: Unable to open the environmental file\n";
        return;
    }
    
    // Prepend model name to title
    ENVFile.List(); ENVFile.Read(Title);
    Title = "bellhopcuda- " + Title;
    
    PRTFile << Title << "\n";
    
    ENVFile.List(); ENVFile.Read(freqinfo->freq0);
    PRTFile << std::setiosflags(std::ios::scientific) << std::setprecision(4);
    PRTFile << " frequency = " << std::setw(11) << freqinfo->freq0 << " Hz\n";
    
    ENVFile.List(); ENVFile.Read(NMedia);
    PRTFile << "Dummy parameter NMedia = " << NMedia << "\n";
    if(NMedia != 1){
        std::cout << "ReadEnvironment: Only one medium or layer is allowed in BELLHOP; sediment layers must be handled using a reflection coefficient\n";
    }
    
    ReadTopOpt(Bdry->Top.hs.Opt, Bdry->Top.hs.bc, FileRoot, 
        ENVFile, PRTFile, ssp, atten);
    
    // *** Top BC ***
    
    if(Bdry->Top.hs.bc == 'A') PRTFile << "   z (m)     alphaR (m/s)   betaR  rho (g/cm^3)  alphaI     betaI\n";
    
    TopBot(freqinfo->freq0, ssp->AttenUnit, fT, Bdry->Top.hs, ENVFile, PRTFile, atten);
    
    // ****** Read in ocean SSP data ******
    
    ENVFile.List(); ENVFile.Read(NPts); ENVFile.Read(Sigma); ENVFile.Read(Bdry->Bot.hs.Depth);
    PRTFile << "\n  Depth = " << std::setw(10) << std::setprecision(2) << Bdry->Bot.hs.Depth << "  m\n";
    
    if(Bdry->Top.hs.Opt[0] == 'A'){
        PRTFile << "Analytic SSP option\n";
        // following is hokey, should be set in Analytic routine
        ssp->NPts = 2;
        ssp->z[0] = RC(0.0);
        ssp->z[1] = Bdry->Bot.hs.Depth;
    }else{
        x = vec2(RC(0.0), Bdry->Bot.hs.Depth);
        InitializeSSP(SSP_CALL_INIT_ARGS);
    }
    
    Bdry->Top.hs.Depth = ssp->z[0]; // Depth of top boundary is taken from first SSP point
    // bottom depth should perhaps be set the same way?
    
    // *** Bottom BC ***
    memcpy(Bdry->Bot.hs.Opt, "  ", 2); // initialize to blanks
    ENVFile.List(); ENVFile.Read(Bdry->Bot.hs.Opt, 6); ENVFile.Read(Sigma);
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
        std::cout << "Unknown bottom option letter in second position\n";
    }
    
    Bdry->Bot.hs.bc = Bdry->Bot.hs.Opt[0];
    TopBot(freqinfo->freq0, ssp->AttenUnit, fT, Bdry->Bot.hs, ENVFile, PRTFile, atten);
    
    // *** source and receiver locations ***
    
    ReadSxSy(false, ENVFile, PRTFile, Pos);
    
    ZMin = Bdry->Top.hs.Depth;
    ZMax = Bdry->Bot.hs.Depth;
    // not sure why I had this
    // ReadSzRz(ZMin + RC(100.0) * (std::nextafter(ZMin, ZMin+RC(1.0)) - ZMin),
    //     ZMax - RC(100.0) * (std::nextafter(ZMax, ZMax+RC(1.0)) - ZMax), 
    //     ENVFile, PRTFile);
    ReadSzRz(ZMin, ZMax, ENVFile, PRTFile, Pos);
    ReadRcvrRanges(ENVFile, PRTFile, Pos);
    ReadfreqVec(Bdry->Top.hs.Opt[5], ENVFile, PRTFile, freqinfo);
    ReadRunType(Beam->RunType, PlotType, ENVFile, PRTFile, Pos);
    
    Depth = ZMax - ZMin; // water depth
    ReadRayElevationAngles(freqinfo->freq0, Depth, Bdry->Top.hs.Opt, Beam->RunType, 
        ENVFile, PRTFile, Angles, Pos);
    
    //LP: Moved to separate function for clarity and modularity.
    ReadBeamInfo(ENVFile, PRTFile, Beam);
} 
