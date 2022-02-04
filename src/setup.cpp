#include "setup.hpp"
#include "readenv.hpp"
#include "output.hpp"

constexpr bool Init_Inline = false;

void setup(std::string FileRoot, 
    std::ofstream &PRTFile, LDOFile &RAYFile, DirectOFile &SHDFile,
    std::string &Title, real &fT,
    BdryType *&Bdry, BdryInfo *&bdinfo, ReflectionInfo *&refl, SSPStructure *&ssp,
    AttenInfo *&atten, Position *&Pos, AnglesStructure *&Angles, FreqInfo *&freqinfo, 
    BeamStructure *&Beam, BeamInfo *&beaminfo, EigenInfo *&eigen, ArrInfo *&arrinfo)
{
    PRTFile.open(FileRoot + ".prt");
    if(!PRTFile.good()){
        std::cout << "Could not open print file: " << FileRoot << ".prt\n";
        std::abort();
    }
    PRTFile << std::unitbuf;
    
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
    eigen = allocate<EigenInfo>();
    arrinfo = allocate<ArrInfo>();
    HSInfo RecycledHS; //Values only initialized once--reused from top to ssp, and ssp to bot
    
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
    eigen->hits = nullptr;
    arrinfo->Arr = nullptr;
    arrinfo->NArr = nullptr;
    
    // Fill in default / "constructor" data
    fT = RL(1.0e20);
    //Bdry: none
    bdinfo->NATIPts = 2;
    bdinfo->NBTYPts = 2;
    memcpy(bdinfo->atiType, "LS", 2);
    memcpy(bdinfo->btyType, "LS", 2);
    //refl: none
    //ssp: none
    atten->t = FL(20.0);
    atten->Salinity = FL(35.0);
    atten->pH = FL(8.0);
    atten->z_bar = FL(0.0);
    Pos->NSx = 1;
    Pos->NSy = 1;
    Angles->Nalpha = 0;
    Angles->Nbeta = 1;
    //LP: not a typo; this is an index, one less than the start of the array,
    //which in Fortran (and in the env file!) is 0. This gets converted to 0-
    //indexed when it is used.
    Angles->iSingle_alpha = 0;
    Angles->iSingle_beta = 0;
    //freqinfo: none
    Beam->epsMultiplier = FL(1.0);
    memcpy(Beam->Type, "G S ", 4);
    //beaminfo: none
    eigen->neigen = 0;
    eigen->memsize = 0;
    arrinfo->MaxNArr = 1;
    RecycledHS.alphaR = FL(1500.0);
    RecycledHS.betaR = FL(0.0);
    RecycledHS.alphaI = FL(0.0);
    RecycledHS.betaI = FL(0.0);
    RecycledHS.rho = FL(1.0);
    
    if(Init_Inline){
        // NPts, Sigma not used by BELLHOP
        Title = PROGRAMNAME "- Calibration case with envfil passed as parameters\n";
        freqinfo->freq0 = FL(250.0);
        // NMedia variable is not used by BELLHOP
        
        // *** Boundary information (type of boundary condition and, if a halfspace, then halfspace info)
        
        memcpy(ssp->AttenUnit, "W", 2); //LP: not a typo--one character string assigned to two
        Bdry->Top.hs.bc    = 'V';
        Bdry->Top.hs.Depth = FL(0.0);
        Bdry->Bot.hs.Depth = FL(100.0);
        memcpy(Bdry->Bot.hs.Opt, "A_", 2);
        Bdry->Bot.hs.bc    = 'A';
        Bdry->Bot.hs.cP    = crci(RL(1.0e20), RL(1590.0), RL(0.5), freqinfo->freq0, freqinfo->freq0,
            ssp->AttenUnit, betaPowerLaw, fT, atten, PRTFile); // compressional wave speed
        Bdry->Bot.hs.cS    = crci(RL(1.0e20), RL(0.0)   , RL(0.0), freqinfo->freq0, freqinfo->freq0,
            ssp->AttenUnit, betaPowerLaw, fT, atten, PRTFile); // shear         wave speed
        Bdry->Bot.hs.rho   = FL(1.2);
        
        // *** sound speed in the water column ***
        
        ssp->Type = 'C'; // interpolation method for SSP
        ssp->NPts = 2;   // number of SSP points
        ssp->z[0]  = FL(0.0);    ssp->z[1]  = FL(100.0);
        ssp->c[0]  = FL(1500.0); ssp->c[1]  = FL(1500.0);
        ssp->cz[0] = FL(0.0);    ssp->cz[1] = FL(0.0); // user should really not have to supply this ...
        
        // *** source and receiver positions ***
        
        Pos->NSz = 1;
        Pos->NRz = 100;
        Pos->NRr = 500;
        
        Pos->Sz = allocate<float>(Pos->NSz); Pos->ws = allocate<float>(Pos->NSz); Pos->iSz = allocate<int32_t>(Pos->NSz);
        Pos->Rz = allocate<float>(Pos->NRz); Pos->wr = allocate<float>(Pos->NRz); Pos->iRz = allocate<int32_t>(Pos->NRz);
        Pos->Rr = allocate<float>(Pos->NRr);
        
        memcpy(Beam->RunType, "C      ", 7);
        memcpy(Beam->Type, "G   ", 4);
        Beam->deltas  = FL(0.0);
        Beam->Box.z   = FL(101.0);
        Beam->Box.r   = FL(5100.0); // meters
        
        Angles->Nalpha = 1789;
        //Angles->alpha = {-80, -70, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80}; // -89 89
        for(int32_t jj=0; jj<Angles->Nalpha; ++jj) Angles->alpha[jj] = (FL(180.0) / Angles->Nalpha) * (real)jj - FL(90.0);
        
        // *** altimetry ***
        
        bdinfo->Top = allocate<BdryPtFull>(2);
        bdinfo->Top[0].x = vec2(-BdryInfinity, RL(0.0));
        bdinfo->Top[1].x = vec2( BdryInfinity, RL(0.0));
        
        ComputeBdryTangentNormal(bdinfo->Top, true, bdinfo);
        
        // *** bathymetry ***
        
        bdinfo->Bot = allocate<BdryPtFull>(2);
        bdinfo->Bot[0].x = vec2(-BdryInfinity, RL(5000.0));
        bdinfo->Bot[1].x = vec2( BdryInfinity, RL(5000.0));
        
        ComputeBdryTangentNormal(bdinfo->Bot, false, bdinfo);
        
        refl->RBot = allocate<ReflectionCoef>(1);
        refl->RTop = allocate<ReflectionCoef>(1);
        //LP: BUG: On this codepath, BELLHOP never sets refl->NBotPts and
        //refl->NTopPts, and they don't have any default value.
        refl->NBotPts = refl->NTopPts = 1;
        
        // *** Source Beam Pattern ***
        beaminfo->NSBPPts = 2;
        beaminfo->SrcBmPat = allocate<real>(2*2);
        beaminfo->SrcBmPat[0*2+0] = FL(-180.0); beaminfo->SrcBmPat[0*2+1] = FL(0.0);
        beaminfo->SrcBmPat[1*2+0] = FL( 180.0); beaminfo->SrcBmPat[1*2+1] = FL(0.0);
        for(int32_t i=0; i<2; ++i) beaminfo->SrcBmPat[i*2+1] = 
            STD::pow(FL(10.0), beaminfo->SrcBmPat[i*2+1] / FL(20.0)); // convert dB to linear scale !!!
    }else{
        ReadEnvironment(FileRoot, PRTFile, Title, fT, Bdry,
            ssp, atten, Pos, Angles, freqinfo, Beam, RecycledHS);
        ReadATI(FileRoot, Bdry->Top.hs.Opt[4], Bdry->Top.hs.Depth, PRTFile, bdinfo); // AlTImetry
        ReadBTY(FileRoot, Bdry->Bot.hs.Opt[1], Bdry->Bot.hs.Depth, PRTFile, bdinfo); // BaThYmetry
        ReadReflectionCoefficient(FileRoot, Bdry->Bot.hs.Opt[0], Bdry->Top.hs.Opt[1], PRTFile, refl); // (top and bottom)
        beaminfo->SBPFlag = Beam->RunType[2];
        ReadPat(FileRoot, PRTFile, beaminfo); // Source Beam Pattern
        // dummy bearing angles
        Pos->Ntheta = 1;
        Pos->theta = allocate<float>(Pos->Ntheta);
        Pos->theta[0] = FL(0.0);
    }
    
    // LP: Moved from WriteHeader
    // receiver bearing angles
    if(Pos->theta == nullptr){
        Pos->theta = allocate<float>(1);
        Pos->theta[0] = FL(0.0); // dummy bearing angle
        Pos->Ntheta = 1;
    }
    // source x-coordinates
    if(Pos->Sx == nullptr){
        Pos->Sx = allocate<float>(1);
        Pos->Sx[0] = FL(0.0); // dummy x-coordinate
        Pos->NSx = 1;
    }
    // source y-coordinates
    if(Pos->Sy == nullptr){
        Pos->Sy = allocate<float>(1);
        Pos->Sy[0] = FL(0.0); // dummy y-coordinate
        Pos->NSy = 1;
    }
    
    OpenOutputFiles(FileRoot, false, Title, Bdry, Pos, Angles, freqinfo, Beam,
        RAYFile, SHDFile);
    
    if(Beam->deltas == FL(0.0)){
        Beam->deltas = (Bdry->Bot.hs.Depth - Bdry->Top.hs.Depth) / FL(10.0); // Automatic step size selection
        PRTFile << "\n Step length,       deltas = " << Beam->deltas << " m (automatically selected)\n";
    }
    
    for(int32_t i=0; i<Angles->Nalpha; ++i) Angles->alpha[i] *= DegRad; // convert to radians
    Angles->Dalpha = FL(0.0);
    if(Angles->Nalpha != 1)
        Angles->Dalpha = (Angles->alpha[Angles->Nalpha-1] - Angles->alpha[0]) / (Angles->Nalpha-1);
    
    // convert range-dependent geoacoustic parameters from user to program units
    if(bdinfo->atiType[1] == 'L'){
        for(int32_t iSeg = 0; iSeg < bdinfo->NATIPts; ++iSeg){
            bdinfo->Top[iSeg].hs.cP = crci(RL(1.0e20), bdinfo->Top[iSeg].hs.alphaR,
                bdinfo->Top[iSeg].hs.alphaI, freqinfo->freq0, freqinfo->freq0,
                {'W', ' '}, betaPowerLaw, fT, atten, PRTFile); // compressional wave speed
            bdinfo->Top[iSeg].hs.cS = crci(RL(1.0e20), bdinfo->Top[iSeg].hs.betaR,
                bdinfo->Top[iSeg].hs.betaI, freqinfo->freq0, freqinfo->freq0,
                {'W', ' '}, betaPowerLaw, fT, atten, PRTFile); // shear         wave speed
        }
    }
    
    if(bdinfo->btyType[1] == 'L'){
        for(int32_t iSeg = 0; iSeg < bdinfo->NBTYPts; ++iSeg){
            bdinfo->Bot[iSeg].hs.cP = crci(RL(1.0e20), bdinfo->Bot[iSeg].hs.alphaR,
                bdinfo->Bot[iSeg].hs.alphaI, freqinfo->freq0, freqinfo->freq0,
                {'W', ' '}, betaPowerLaw, fT, atten, PRTFile); // compressional wave speed
            bdinfo->Bot[iSeg].hs.cS = crci(RL(1.0e20), bdinfo->Bot[iSeg].hs.betaR,
                bdinfo->Bot[iSeg].hs.betaI, freqinfo->freq0, freqinfo->freq0,
                {'W', ' '}, betaPowerLaw, fT, atten, PRTFile); // shear         wave speed
        }
    }
    
    PRTFile << "\n";
}
