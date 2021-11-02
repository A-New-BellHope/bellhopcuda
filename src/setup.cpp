#include "setup.hpp"
#include "readenv.hpp"
#include "output.hpp"

constexpr bool Init_Inline = false;

void setup(int argc, char **argv, 
    std::ofstream &PRTFile, std::ofstream &RAYFile, std::ofstream &ARRFile,
    std::string &Title, real &fT,
    BdryType *&Bdry, BdryInfo *&bdinfo, ReflectionInfo *&refl, SSPStructure *&ssp,
    AttenInfo *&atten, Position *&Pos, AnglesStructure *&Angles, FreqInfo *&freqinfo, 
    BeamStructure *&Beam, BeamInfo *&beaminfo)
{
    // Command-line parameters
    if(argc < 2){
        std::cout << "Must provide FileRoot as command-line parameter\n";
        std::abort();
    }
    std::string FileRoot(argv[1]);
    PRTFile.open(WithExtension(FileRoot, ".prt"));
    if(!PRTFile.good()){
        std::cout << "Could not open print file: " << WithExtension(FileRoot, ".prt") << "\n";
        std::abort();
    }
    
    // Allocate main structs
    Bdry = allocate<BdryType>();
    bdinfo = allocate<BdryInfo>();
    refl = allocate<ReflectionInfo>();
    ssp = allocate<SSPStructure>();
    atten = allocate<AttenInfo>();
    Pos = allocate<Position>();
    Angles = allocate<AnglesStructure>();
    freqinfo = allocate<FreqInfo>();
    Beam = allocate<BeamStructure>();
    beaminfo = allocate<BeamInfo>();
    
    // Debugging: Fill structs with garbage data to help detect uninitialized vars
    memset(Bdry, 0xFE, sizeof(BdryType));
    memset(bdinfo, 0xFE, sizeof(BdryInfo));
    memset(refl, 0xFE, sizeof(ReflectionInfo));
    memset(ssp, 0xFE, sizeof(SSPStructure));
    memset(atten, 0xFE, sizeof(AttenInfo));
    memset(Pos, 0xFE, sizeof(AttenInfo));
    memset(Angles, 0xFE, sizeof(AnglesStructure));
    memset(freqinfo, 0xFE, sizeof(FreqInfo));
    memset(Beam, 0xFE, sizeof(BeamStructure));
    memset(beaminfo, 0xFE, sizeof(BeamInfo));
    
    // Set pointers to null because BELLHOP checks if some of them are allocated
    // before allocating them
    bdinfo->Top = nullptr;
    bdinfo->Bot = nullptr;
    refl->RBot = nullptr;
    refl->RTop = nullptr;
    ssp->cMat = nullptr;
    ssp->czMat = nullptr;
    ssp->cMat3 = nullptr;
    ssp->czMat3 = nullptr;
    Pos->iSz = nullptr;
    Pos->iRz = nullptr;
    Pos->Sx = nullptr;
    Pos->Sy = nullptr;
    Pos->Sz = nullptr;
    Pos->Rr = nullptr;
    Pos->Rz = nullptr;
    Pos->ws = nullptr;
    Pos->wr = nullptr;
    Pos->theta = nullptr;
    Angles->alpha = nullptr;
    Angles->beta = nullptr;
    freqinfo->freqVec = nullptr;
    beaminfo->SrcBmPat = nullptr;
    
    // Fill in default / "constructor" data
    fT = RC(1.0e20);
    //Bdry: none
    bdinfo->NATIPts = 2;
    bdinfo->NBTYPts = 2;
    memcpy(bdinfo->atiType, "LS", 2);
    memcpy(bdinfo->btyType, "LS", 2);
    //refl: none
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
    //freqinfo: none
    Beam->epsMultiplier = RC(1.0);
    memcpy(Beam->Type, "G S ", 4);
    //beaminfo: none
    
    if(Init_Inline){
        // NPts, Sigma not used by BELLHOP
        Title = "bellhopcuda- Calibration case with envfil passed as parameters\n";
        freqinfo->freq0 = RC(250.0);
        // NMedia variable is not used by BELLHOP
        
        // *** Boundary information (type of boundary condition and, if a halfspace, then halfspace info)
        
        memcpy(ssp->AttenUnit, "W", 2); //LP: not a typo--one character string assigned to two
        Bdry->Top.hs.bc    = 'V';
        Bdry->Top.hs.Depth = RC(0.0);
        Bdry->Bot.hs.Depth = RC(100.0);
        memcpy(Bdry->Bot.hs.Opt, "A_", 2);
        Bdry->Bot.hs.bc    = 'A';
        Bdry->Bot.hs.cP    = crci(RC(1.0e20), RC(1590.0), RC(0.5), freqinfo->freq0, freqinfo->freq0,
            ssp->AttenUnit, betaPowerLaw, fT, atten, PRTFile); // compressional wave speed
        Bdry->Bot.hs.cS    = crci(RC(1.0e20), RC(0.0)   , RC(0.0), freqinfo->freq0, freqinfo->freq0,
            ssp->AttenUnit, betaPowerLaw, fT, atten, PRTFile); // shear         wave speed
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
        
        Pos->Sz = allocate<real>(Pos->NSz); Pos->ws = allocate<real>(Pos->NSz); Pos->iSz = allocate<int32_t>(Pos->NSz);
        Pos->Rz = allocate<real>(Pos->NRz); Pos->wr = allocate<real>(Pos->NRz); Pos->iRz = allocate<int32_t>(Pos->NRz);
        Pos->Rr = allocate<real>(Pos->NRr);
        
        memcpy(Beam->RunType, "C      ", 7);
        memcpy(Beam->Type, "G   ", 4);
        Beam->deltas  = RC(0.0);
        Beam->Box.z   = RC(101.0);
        Beam->Box.r   = RC(5100.0); // meters
        
        Angles->Nalpha = 1789;
        //Angles->alpha = {-80, -70, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80}; // -89 89
        for(int32_t jj=0; jj<Angles->Nalpha; ++jj) Angles->alpha[jj] = (RC(180.0) / Angles->Nalpha) * (real)jj - RC(90.0);
        
        // *** altimetry ***
        
        bdinfo->Top = allocate<BdryPtFull>(2);
        bdinfo->Top[0].x = vec2(-BdryInfinity(), RC(0.0));
        bdinfo->Top[1].x = vec2( BdryInfinity(), RC(0.0));
        
        ComputeBdryTangentNormal(bdinfo->Top, true, bdinfo);
        
        // *** bathymetry ***
        
        bdinfo->Bot = allocate<BdryPtFull>(2);
        bdinfo->Bot[0].x = vec2(-BdryInfinity(), RC(5000.0));
        bdinfo->Bot[1].x = vec2( BdryInfinity(), RC(5000.0));
        
        ComputeBdryTangentNormal(bdinfo->Bot, false, bdinfo);
        
        refl->RBot = allocate<ReflectionCoef>(1);
        refl->RTop = allocate<ReflectionCoef>(1);
        //LP: BUG: On this codepath, BELLHOP never sets refl->NBotPts and
        //refl->NTopPts, and they don't have any default value.
        refl->NBotPts = refl->NTopPts = 1;
        
        // *** Source Beam Pattern ***
        beaminfo->NSBPPts = 2;
        beaminfo->SrcBmPat = allocate<real>(2*2);
        beaminfo->SrcBmPat[0*2+0] = RC(-180.0); beaminfo->SrcBmPat[0*2+1] = RC(0.0);
        beaminfo->SrcBmPat[1*2+0] = RC( 180.0); beaminfo->SrcBmPat[1*2+1] = RC(0.0);
        for(int32_t i=0; i<2; ++i) beaminfo->SrcBmPat[i*2+1] = 
            STD::pow(RC(10.0), beaminfo->SrcBmPat[i*2+1] / RC(20.0)); // convert dB to linear scale !!!
    }else{
        ReadEnvironment(FileRoot, PRTFile, Title, fT, Bdry,
            ssp, atten, Pos, Angles, freqinfo, Beam);
        ReadATI(FileRoot, Bdry->Top.hs.Opt[4], Bdry->Top.hs.Depth, PRTFile, bdinfo); // AlTImetry
        ReadBTY(FileRoot, Bdry->Bot.hs.Opt[1], Bdry->Bot.hs.Depth, PRTFile, bdinfo); // BaThYmetry
        ReadReflectionCoefficient(FileRoot, Bdry->Bot.hs.Opt[0], Bdry->Top.hs.Opt[1], PRTFile, refl); // (top and bottom)
        beaminfo->SBPFlag = Beam->RunType[2];
        ReadPat(FileRoot, PRTFile, beaminfo); // Source Beam Pattern
        // dummy bearing angles
        Pos->Ntheta = 1;
        Pos->theta = allocate<real>(Pos->Ntheta);
        Pos->theta[0] = RC(0.0);
    }
    
    OpenOutputFiles(FileRoot, false, Title, Bdry, Pos, Angles, freqinfo, Beam,
        RAYFile, ARRFile);
}

void core_setup(std::ofstream &PRTFile, const real &fT,
    const BdryType *Bdry, const BdryInfo *bdinfo, const AttenInfo *atten, 
    AnglesStructure *Angles, const FreqInfo *freqinfo, BeamStructure *Beam/*, 
    InfluenceInfo *inflinfo, ArrivalsInfo *arrinfo*/)
{
    if(Beam->deltas == RC(0.0)){
        Beam->deltas = (Bdry->Bot.hs.Depth - Bdry->Top.hs.Depth) / RC(10.0); // Automatic step size selection
        PRTFile << "\n Step length,       deltas = " << Beam->deltas << " m (automatically selected)\n";
    }
    
    for(int32_t i=0; i<Angles->Nalpha; ++i) Angles->alpha[i] *= DegRad; // convert to radians
    Angles->Dalpha = RC(0.0);
    if(Angles->Nalpha != 1)
        Angles->Dalpha = (Angles->alpha[Angles->Nalpha-1] - Angles->alpha[0]) / (Angles->Nalpha-1);
    
    // convert range-dependent geoacoustic parameters from user to program units
    if(bdinfo->atiType[1] == 'L'){
        for(int32_t iSeg = 0; iSeg < bdinfo->NATIPts; ++iSeg){
            bdinfo->Top[iSeg].hs.cP = crci(RC(1.0e20), bdinfo->Top[iSeg].hs.alphaR,
                bdinfo->Top[iSeg].hs.alphaI, freqinfo->freq0, freqinfo->freq0,
                {'W', ' '}, betaPowerLaw, fT, atten, PRTFile); // compressional wave speed
            bdinfo->Top[iSeg].hs.cS = crci(RC(1.0e20), bdinfo->Top[iSeg].hs.betaR,
                bdinfo->Top[iSeg].hs.betaI, freqinfo->freq0, freqinfo->freq0,
                {'W', ' '}, betaPowerLaw, fT, atten, PRTFile); // shear         wave speed
        }
    }
    
    if(bdinfo->btyType[1] == 'L'){
        for(int32_t iSeg = 0; iSeg < bdinfo->NBTYPts; ++iSeg){
            bdinfo->Bot[iSeg].hs.cP = crci(RC(1.0e20), bdinfo->Bot[iSeg].hs.alphaR,
                bdinfo->Bot[iSeg].hs.alphaI, freqinfo->freq0, freqinfo->freq0,
                {'W', ' '}, betaPowerLaw, fT, atten, PRTFile); // compressional wave speed
            bdinfo->Bot[iSeg].hs.cS = crci(RC(1.0e20), bdinfo->Bot[iSeg].hs.betaR,
                bdinfo->Bot[iSeg].hs.betaI, freqinfo->freq0, freqinfo->freq0,
                {'W', ' '}, betaPowerLaw, fT, atten, PRTFile); // shear         wave speed
        }
    }
    
    /*
    if(Beam->RunType[4] == 'I'){
        inflinfo->NRz_per_range = 1; // irregular grid
    }else{
        inflinfo->NRz_per_range = Pos->NRz; // rectilinear grid
    }
    
    // for a TL calculation, allocate space for the pressure matrix
    // TODO need a separate one for each source coord?
    switch(Beam->RunType[0]){
    case 'C':
    case 'S':
    case 'I':
        // TL calculation
        inflinfo->u = allocate<real>(inflinfo->NRz_per_range * Pos->NRr);
        break;
    case 'A':
    case 'a':
    case 'R':
    case 'E':
        // Arrivals calculation
        inflinfo->u = allocate<real>(1); // open a dummy variable
        break;
    default:
        std::cout << "RunType[0] == " << Beam->RunType[0] << " case not handled by BELLHOP\n";
        std::abort();
    }
    
    switch(Beam->RunType[0]){
    case 'A':
    case 'a':
        arrinfo->MaxNArr = std::max(ArrivalsStorage / (inflinfo->NRz_per_range * Pos->NRr), MinNArr);
        PRTFile << "\n( Maximum # of arrivals = " << arrinfo->MaxNArr << ")\n";
        break;
    default:
        arrinfo->MaxNArr = 1;
        arrinfo->arr = allocate<TODO>(inflinfo->NRz_per_range * Pos->NRr * 1);
        arrinfo->NArr = allocate<int32_t>(inflinfo->NRz_per_range * Pos->NRr);
    }
    
    memset(arrinfo->NArr, 0, inflinfo->NRz_per_range * Pos->NRr * sizeof(int32_t));
    */
    
    PRTFile << "\n";
}
