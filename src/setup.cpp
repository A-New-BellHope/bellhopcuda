#include "common.hpp"

constexpr bool Init_Inline = false;

void setup(int argc, char **argv,
    FreqInfo *&freqinfo, BdryType *&Bdry, BdryInfo *&binfo, SSPStructure *&ssp,
    AttenInfo *&atten, Position *&Pos, AnglesStructure *&Angles, BeamStructure *&Beam
    )
{
    // Command-line parameters
    if(argc < 2){
        std::cout << "Must provide FileRoot as command-line parameter\n";
        std::abort();
    }
    std::string FileRoot(argv[1]);
    std::ostream PRTFile(WithExtension(FileRoot, ".prt"));
    if(!PRTFile.good()){
        std::cout << "Could not open print file: " << WithExtension(FileRoot, ".prt") << "\n";
        std::abort();
    }
    std::string Title;
    real fT;
    
    // Allocate main structs
    freqinfo = allocate<FreqInfo>();
    Bdry = allocate<BdryType>();
    binfo = allocate<BdryInfo>();
    ssp = allocate<SSPStructure>();
    atten = allocate<AttenInfo>();
    Pos = allocate<Position>();
    Angles = allocate<AnglesStructure>();
    Beam = allocate<BeamStructure>();
    
    // Debugging: Fill structs with garbage data to help detect uninitialized vars
    memset(freqinfo, 0xFE, sizeof(FreqInfo));
    memset(Bdry, 0xFE, sizeof(BdryType));
    memset(binfo, 0xFE, sizeof(BdryInfo));
    memset(ssp, 0xFE, sizeof(SSPStructure));
    memset(atten, 0xFE, sizeof(AttenInfo));
    memset(Pos, 0xFE, sizeof(AttenInfo));
    memset(Angles, 0xFE, sizeof(AnglesStructure));
    memset(Beam, 0xFE, sizeof(BeamStructure));
    
    //Fill in default / "constructor" data
    fT = RC(1.0e20);
    //freqinfo: none
    //Bdry: none
    binfo->NATIPts = 2;
    binfo->NBTYPts = 2;
    binfo->atiType = {'L', 'S'};
    binfo->btyType = {'L', 'S'};
    //ssp: none
    atten->t = RC(20.0);
    atten->Salinity = RC(35.0);
    atten->pH = RC(8.0);
    atten->z_bar = RC(0.0);
    Pos->NSx = 1;
    Pos->NSy = 1;
    Angles->Nalpha = 0;
    Angles->Nbeta = 1;
    Angles->iSingle_alpha = -1; //LP: not a typo; this is an index, one less than the
    Angles->iSingle_beta = -1; //start of the array, which in Fortran is 0 but C++ is -1
    Beam->epsMultiplier = RC(1.0);
    Beam->Type = {'G', ' ', 'S', ' '};
    
    if(Init_Inline){
        // NPts, Sigma not used by BELLHOP
        Title = "bellhopcuda- Calibration case with envfil passed as parameters\n";
        freqinfo->freq0 = RC(250.0);
        // NMedia variable is not used by BELLHOP
        
        // *** Boundary information (type of boundary condition and, if a halfspace, then halfspace info)
        
        ssp->AttenUnit     = {'W', '\0'};
        Bdry->Top.hs.bc    = 'V';
        Bdry->Top.hs.Depth = RC(0.0);
        Bdry->Bot.hs.Depth = RC(100.0);
        Bdry->Bot.hs.Opt   = {'A', '_'};
        Bdry->Bot.hs.bc    = 'A';
        Bdry->Bot.hs.cp    = crci(RC(1.0e20), RC(1590.0), RC(0.5), freqinfo->freq0, freqinfo->freq0,
            ssp->AttenUnit, betaPowerLaw, fT); // compressional wave speed
        Bdry->Bot.hs.cs    = crci(RC(1.0e20), RC(0.0)   , RC(0.0), freqinfo->freq0, freqinfo->freq0,
            ssp->AttenUnit, betaPowerLaw, fT); // shear         wave speed
        Bdry->Bot.hs.rho   = RC(1.2);
        
        // *** sound speed in the water column ***
        
        ssp->Type = 'C'; // interpolation method for SSP
        ssp->NPts = 2;   // number of SSP points
        ssp->z[0]  = RC(0.0);    ssp->z[1]  = RC(100.0);
        ssp->c[0]  = RC(1500.0); ssp->c[1]  = RC(1500.0);
        ssp->cz[0] = RC(0.0);    ssp->cz[1] = RC(0.0); // user should really not have to supply this ...
        
        // *** source and receiver positions ***
        
        Pos->NSz = 1;
        Pos->NRz = 100;
        Pos->NRr = 500;
        
        Pos->sz = allocate<real>(Pos->NSz); Pos->ws = allocate<real>(Pos->NSz); Pos->isz = allocate<int32_t>(Pos->NSz);
        Pos->rz = allocate<real>(Pos->NRz); Pos->wr = allocate<real>(Pos->NRz); Pos->irz = allocate<int32_t>(Pos->NRz);
        Pos->Rr = allocate<real>(Pos->NRr);
        
        Beam->RunType = 'C';
        Beam->Type    = {'G', ' ', ' ', ' '};
        Beam->deltas  = RC(0.0);
        Beam->Box.z   = RC(101.0);
        Beam->Box.r   = RC(5100.0); // meters
        
        Angles->Nalpha = 1789;
        //Angles->alpha = {-80, -70, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80}; // -89 89
        for(int32_t jj=0; jj<Angles->Nalpha; ++jj) Angles->alpha[jj] = (RC(180.0) / Angles->Nalpha) * (real)jj - RC(90.0);
        
        // *** altimetry ***
        
        binfo->Top = allocate<BdryPtFull>(2);
        binfo->Top[0].x = vec2(-BdryInfinity(), RC(0.0));
        binfo->Top[1].x = vec2( BdryInfinity(), RC(0.0));
        
        ComputeBdryTangentNormal(binfo->Top, true, binfo);
        
        // *** bathymetry ***
        
        binfo->Bot = allocate<BdryPtFull>(2);
        binfo->Bot[0].x = vec2(-BdryInfinity(), RC(5000.0));
        binfo->Bot[1].x = vec2( BdryInfinity(), RC(5000.0));
        
        ComputeBdryTangentNormal(binfo->Bot, false, binfo);
        
        //TODO
    }
}
