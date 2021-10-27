#include "common.hpp"

constexpr bool Init_Inline = false;

void setup(int argc, char **argv,
    FreqInfo *&freqinfo, BdryType *&Bdry, SSPStructure *&ssp, AttenInfo *&atten, 
    Position *&Pos, AnglesStructure *&Angles, BeamStructure *&Beam
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
    real fT = RC(1.0e20);
    
    // Allocate main structs
    freqinfo = allocate<FreqInfo>();
    Bdry = allocate<BdryType>();
    ssp = allocate<SSPStructure>();
    atten = allocate<AttenInfo>();
    Pos = allocate<Position>();
    Angles = allocate<AnglesStructure>();
    Beam = allocate<BeamStructure>();
    
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
        
        //TODO
    }
}
