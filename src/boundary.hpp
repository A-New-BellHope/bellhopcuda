#pragma once
#include "common.hpp"
#include "ldio.hpp"
#include "ssp.hpp"

struct BdryPtFull {
    vec2 x, t, n; // coordinate, tangent, and outward normal for a segment
    vec2 Nodet, Noden; // tangent and normal at the node, if the curvilinear option is used
    real Len, kappa; // length and curvature of a segement
    real Dx, Dxx, Dss; // first, second derivatives wrt depth; s is along tangent
    HSInfo hs;
};
 
struct BdryInfo {
    int32_t NATIPts, NBTYPts;
    char atiType[2];
    char btyType[2];
    BdryPtFull *Top, *Bot;
};

constexpr int32_t Bdry_Number_to_Echo = 21;

constexpr HOST_DEVICE inline real BdryInfinity(){
    return STD::sqrt(REAL_MAX) / RC(1.0e5);
}

/**
 * Get the Top segment info (index and range interval) for range, r
 *
 * rTopSeg: segment limits in range
 */
HOST_DEVICE inline void GetTopSeg(real r, int32_t &IsegTop, vec2 &rTopSeg,
    const BdryInfo *bdinfo)
{
    // LP: bdinfo->Top.x is checked for being monotonic at load time, so we can
    // linearly search out from the last position, usually only have to move
    // by 1
    int32_t n = bdinfo->NATIPts;
    IsegTop = STD::max(IsegTop, 0);
    IsegTop = STD::min(IsegTop, n-2);
    while(IsegTop >= 0 && bdinfo->Top[IsegTop].x.x > r) --IsegTop;
    while(IsegTop >= 0 && IsegTop < n-1 && bdinfo->Top[IsegTop+1].x.x < r) ++IsegTop;
    if(IsegTop < 0 || IsegTop >= n-1){
        // IsegTop MUST LIE IN [0, NATIPts-2]
        printf("Error: GetTopSeg: Top altimetry undefined above the ray, r=%g\n", r);
        bail();
    }
    rTopSeg.x = bdinfo->Top[IsegTop].x.x;
    rTopSeg.y = bdinfo->Top[IsegTop+1].x.x;
}

/**
 * Get the Bottom segment info (index and range interval) for range, r
 *
 * rBotSeg: segment limits in range
 */
HOST_DEVICE inline void GetBotSeg(real r, int32_t &IsegBot, vec2 &rBotSeg,
    const BdryInfo *bdinfo)
{
    // LP: bdinfo->Bot.x is checked for being monotonic at load time, so we can
    // linearly search out from the last position, usually only have to move
    // by 1
    int32_t n = bdinfo->NBTYPts;
    IsegBot = STD::max(IsegBot, 0);
    IsegBot = STD::min(IsegBot, n-2);
    while(IsegBot >= 0 && bdinfo->Bot[IsegBot].x.x > r) --IsegBot;
    while(IsegBot >= 0 && IsegBot < n-1 && bdinfo->Bot[IsegBot+1].x.x < r) ++IsegBot;
    if(IsegBot < 0 || IsegBot >= n-1){
        // IsegBot MUST LIE IN [0, NBTYPts-2]
        printf("Error: GetBotSeg: Bottom bathymetry undefined below the source, r=%g\n", r);
        bail();
    }
    rBotSeg.x = bdinfo->Bot[IsegBot].x.x;
    rBotSeg.y = bdinfo->Bot[IsegBot+1].x.x;
}

/**
 * Does some pre-processing on the boundary points to pre-compute segment
 * lengths  (%Len),
 * tangents (%t, %nodet),
 * normals  (%n, %noden), and
 * curvatures (%kappa)
 *
 * The boundary is also extended with a constant depth to infinity to cover cases where the ray
 * exits the domain defined by the user
 */
inline void ComputeBdryTangentNormal(BdryPtFull *Bdry, bool isTop, BdryInfo *bdryInfo)
{
    int32_t NPts = 0;
    real *phi;
    real sss;
    char CurvilinearFlag = '-';
    
    if(!isTop){
        NPts = bdryInfo->NBTYPts;
        CurvilinearFlag = bdryInfo->btyType[0];
    }else{
        NPts = bdryInfo->NATIPts;
        CurvilinearFlag = bdryInfo->atiType[0];
    }
    
    // extend the bathymetry to +/- infinity in a piecewise constant fashion

    Bdry[0     ].x[0] = -BdryInfinity();
    Bdry[0     ].x[1] = Bdry[1     ].x[1];
    Bdry[0     ].hs   = Bdry[1     ].hs;
    Bdry[NPts-1].x[0] =  BdryInfinity();
    Bdry[NPts-1].x[1] = Bdry[NPts-2].x[1];
    Bdry[NPts-1].hs   = Bdry[NPts-2].hs;
    
    // compute tangent and outward-pointing normal to each bottom segment
    // tBdry[0][:] = xBdry[0][1:NPts-1] - xBdry[0][0:NPts-2]
    // tBdry[1][:] = xBdry[1][1:NPts-1] - xBdry[1][0:NPts-2]
    // above caused compiler problems
    // LP: C++ obviously does not have vector slicing, but you get the idea.
    
    for(int32_t ii=0; ii<NPts-1; ++ii){
        Bdry[ii].t  = Bdry[ii+1].x - Bdry[ii].x;
        Bdry[ii].Dx = Bdry[ii].t[1] / Bdry[ii].t[0]; // first derivative
        // std::cout << "Dx, t " << Bdry[ii].Dx << " " << Bdry[ii].x << " " << (RC(1.0) / (Bdry[ii].x[1] / RC(500.0)) << "\n";
        
        // normalize the tangent vector
        Bdry[ii].Len = glm::length(Bdry[ii].t);
        Bdry[ii].t  *= RC(1.0) / Bdry[ii].Len;
        
        Bdry[ii].n[0] = (isTop ? RC(1.0) : RC(-1.0)) * Bdry[ii].t[1];
        Bdry[ii].n[1] = (isTop ? RC(-1.0) : RC(1.0)) * Bdry[ii].t[0];
    }
    
    if(CurvilinearFlag == 'C'){
        // curvilinear option: compute tangent and normal at node by averaging normals on adjacent segments
        // averaging two centered differences is equivalent to forming a single centered difference of two steps ...
        for(int32_t ii=1; ii<NPts-1; ++ii){
            sss = Bdry[ii-1].Len / (Bdry[ii-1].Len + Bdry[ii].Len);
            sss = RC(0.5); //LP: ???
            Bdry[ii].Nodet = (RC(1.0) - sss) * Bdry[ii-1].t + sss * Bdry[ii].t;
        }
        
        Bdry[0     ].Nodet = vec2(RC(1.0), RC(0.0)); // tangent left-end  node
        Bdry[NPts-1].Nodet = vec2(RC(1.0), RC(0.0)); // tangent right-end node
        
        for(int32_t ii=0; ii<NPts; ++ii){
            Bdry[ii].Noden[0] = (isTop ? RC(1.0) : RC(-1.0)) * Bdry[ii].Nodet[1];
            Bdry[ii].Noden[1] = (isTop ? RC(-1.0) : RC(1.0)) * Bdry[ii].Nodet[0];
        }
        
        // compute curvature in each segment
        phi = allocate<real>(NPts);
        // this is the angle at each node
        for(int32_t i=0; i<NPts; ++i) phi[i] = STD::atan2(Bdry[i].Nodet[1], Bdry[i].Nodet[0]);
        
        for(int32_t ii=0; ii<NPts-1; ++ii){
            Bdry[ii].kappa = (phi[ii+1] - phi[ii]) / Bdry[ii].Len; // this is curvature = dphi/ds
            Bdry[ii].Dxx   = (Bdry[ii+1].Dx - Bdry[ii].Dx) / // second derivative
                             (Bdry[ii+1].x[0] - Bdry[ii].x[0]); 
            Bdry[ii].Dss   = Bdry[ii].Dxx * CUBE(Bdry[ii].t[0]); // derivative in direction of tangent
            //std::cout << "kappa, Dss, Dxx " << Bdry[ii].kappa << " " << Bdry[ii].Dss << " " << Bdry[ii].Dxx
            //    << " " << RC(1.0) / ((RC(8.0) / SQ(RC(1000.0))) * CUBE(STD::abs(Bdry[ii].x[1])))
            //    << " " << Bdry[ii].x[1] << " "
            //    << RC(-1.0) / (RC(4.0) * CUBE(Bdry[ii].x[1]) / RC(1000000.0))
            //    << " " << Bdry[ii].x[1] << "\n";
            
            Bdry[ii].kappa = Bdry[ii].Dss; // over-ride kappa !!!!!
        }
    }else{
        for(int32_t i=0; i<NPts; ++i) Bdry[i].kappa = RC(0.0);
    }
}

inline void ReadATI(std::string FileRoot, char TopATI, real DepthT,
    std::ofstream &PRTFile, BdryInfo *bdinfo)
{
    //LP: Removed phi, which got allocated but not used. Just to check that
    //there was enough memory for it to be allocated later? Or bug?
    switch(TopATI){
    case '~':
    case '*':{
        PRTFile << "__________________________________________________________________________\n\n";
        PRTFile << "Using top-altimetry file\n";
        
        LDIFile ATIFile(FileRoot + ".ati");
        if(!ATIFile.Good()){
            PRTFile << "ATIFile = " << FileRoot << ".ati\n";
            std::cout << "ReadATI: Unable to open altimetry file\n";
            std::abort();
        }
        
        LIST(ATIFile); ATIFile.Read(bdinfo->atiType, 2);
        switch(bdinfo->atiType[0]){
        case 'C':
            PRTFile << "Curvilinear Interpolation\n"; break;
        case 'L':
            PRTFile << "Piecewise linear interpolation\n"; break;
        default:
            std::cout << "ReadATI: Unknown option for selecting altimetry interpolation\n";
            std::abort();
        }
        
        LIST(ATIFile); ATIFile.Read(bdinfo->NATIPts);
        PRTFile << "Number of altimetry points = " << bdinfo->NATIPts << "\n";
        bdinfo->NATIPts += 2; // we'll be extending the altimetry to infinity to the left and right
        
        bdinfo->Top = allocate<BdryPtFull>(bdinfo->NATIPts);
        
        PRTFile << "\n Range (km)  Depth (m)\n";
        
        for(int32_t ii=1; ii<bdinfo->NATIPts-1; ++ii){
            switch(bdinfo->atiType[1]){
            case 'S':
            case ' ':
                LIST(ATIFile); ATIFile.Read(bdinfo->Top[ii].x);
                if(ii < Bdry_Number_to_Echo || ii == bdinfo->NATIPts){ // echo some values
                    //LP: Bug? Even in the Fortran, ii should never equal
                    //NATIPts as the loop ends at NATIPts-1.
                    PRTFile << std::setprecision(3) << bdinfo->Top[ii].x << "\n";
                }
                break;
            case 'L':
                LIST(ATIFile); ATIFile.Read(bdinfo->Top[ii].x);
                ATIFile.Read(bdinfo->Top[ii].hs.alphaR);
                ATIFile.Read(bdinfo->Top[ii].hs.betaR);
                ATIFile.Read(bdinfo->Top[ii].hs.rho);
                ATIFile.Read(bdinfo->Top[ii].hs.alphaI);
                ATIFile.Read(bdinfo->Top[ii].hs.betaI);
                if(ii < Bdry_Number_to_Echo || ii == bdinfo->NATIPts){ // echo some values
                    PRTFile << std::setprecision(3) << bdinfo->Top[ii].x << " "
                        << bdinfo->Top[ii].hs.alphaR << " "
                        << bdinfo->Top[ii].hs.betaR << " "
                        << bdinfo->Top[ii].hs.rho << " "
                        << bdinfo->Top[ii].hs.alphaI << " "
                        << bdinfo->Top[ii].hs.betaI << "\n";
                }
                break;
            default:
                std::cout << "ReadATI: Unknown option for selecting altimetry option\n";
                std::abort();
            }
            
            if(bdinfo->Top[ii].x[1] < DepthT){
                std::cout << "BELLHOP:ReadATI: Altimetry rises above highest point in the sound speed profile\n";
                std::abort();
            }
        }
        
        // Convert ranges in km to m
        for(int32_t i=0; i<bdinfo->NATIPts; ++i) bdinfo->Top[i].x[0] *= RC(1000.0);
        
        }break;
    default:
        bdinfo->Top = allocate<BdryPtFull>(2);
        bdinfo->Top[0].x = vec2(-STD::sqrt(REAL_MAX) / RC(1.0e5), DepthT);
        bdinfo->Top[1].x = vec2( STD::sqrt(REAL_MAX) / RC(1.0e5), DepthT);
    }
    
    ComputeBdryTangentNormal(bdinfo->Top, true, bdinfo);
    
    if(!monotonic(&bdinfo->Top[0].x.x, bdinfo->NATIPts, sizeof(BdryPtFull)/sizeof(real), 0)){
        std::cout << "BELLHOP:ReadATI: Altimetry ranges are not monotonically increasing\n";
        std::abort();
    }
}

inline void ReadBTY(std::string FileRoot, char BotBTY, real DepthB,
    std::ofstream &PRTFile, BdryInfo *bdinfo)
{
    switch(BotBTY){
    case '~':
    case '*':{
        PRTFile << "__________________________________________________________________________\n\n";
        PRTFile << "Using bottom-bathymetry file\n";
        
        LDIFile BTYFile(FileRoot + ".bty");
        if(!BTYFile.Good()){
            PRTFile << "BTYFile = " << FileRoot << ".bty\n";
            std::cout << "ReadATI: Unable to open bathymetry file\n";
            std::abort();
        }
        
        LIST(BTYFile); BTYFile.Read(bdinfo->btyType, 2);
        
        switch(bdinfo->btyType[0]){
        case 'C':
            PRTFile << "Curvilinear Interpolation\n"; break;
        case 'L':
            PRTFile << "Piecewise linear interpolation\n"; break;
        default:
            std::cout << "ReadBTY: Unknown option for selecting bathymetry interpolation\n";
            std::abort();
        }
        
        LIST(BTYFile); BTYFile.Read(bdinfo->NBTYPts);
        PRTFile << "Number of bathymetry points = " << bdinfo->NBTYPts << "\n";
        
        bdinfo->NBTYPts += 2; // we'll be extending the bathymetry to infinity to the left and right
        bdinfo->Bot = allocate<BdryPtFull>(bdinfo->NBTYPts);
        
        PRTFile << "\n";
        switch(bdinfo->btyType[1]){
        case 'S':
        case ' ':
            PRTFile << "Short format (bathymetry only)\n";
            PRTFile << " Range (km)  Depth (m)\n";
            break;
        case 'L':
            PRTFile << "Long format (bathymetry and geoacoustics)\n";
            PRTFile << "Range (km)  Depth (m)  alphaR (m/s)  betaR  rho (g/cm^3)  alphaI     betaI\n";
            break;
        default:
            std::cout << "ReadBTY: Unknown option for selecting bathymetry interpolation\n";
            std::abort();
        }
        
        for(int32_t ii=1; ii<bdinfo->NBTYPts-1; ++ii){
            switch(bdinfo->btyType[1]){
            case 'S':
            case ' ':
                LIST(BTYFile); BTYFile.Read(bdinfo->Bot[ii].x);
                if(ii < Bdry_Number_to_Echo || ii == bdinfo->NBTYPts){ // echo some values
                    //LP: Bug? Even in the Fortran, ii should never equal
                    //NBTYPts as the loop ends at NBTYPts-1.
                    PRTFile << std::setprecision(3) << bdinfo->Bot[ii].x << "\n";
                }
                break;
            case 'L':
                LIST(BTYFile); BTYFile.Read(bdinfo->Bot[ii].x);
                BTYFile.Read(bdinfo->Bot[ii].hs.alphaR);
                BTYFile.Read(bdinfo->Bot[ii].hs.betaR);
                BTYFile.Read(bdinfo->Bot[ii].hs.rho);
                BTYFile.Read(bdinfo->Bot[ii].hs.alphaI);
                BTYFile.Read(bdinfo->Bot[ii].hs.betaI);
                if(ii < Bdry_Number_to_Echo || ii == bdinfo->NBTYPts){ // echo some values
                    PRTFile << std::setprecision(2) << bdinfo->Bot[ii].x << " "
                        << bdinfo->Bot[ii].hs.alphaR << " "
                        << bdinfo->Bot[ii].hs.betaR << " "
                        << bdinfo->Bot[ii].hs.rho << " " << std::setprecision(4)
                        << bdinfo->Bot[ii].hs.alphaI << " "
                        << bdinfo->Bot[ii].hs.betaI << "\n";
                }
                break;
            default:
                std::cout << "ReadBTY: Unknown option for selecting bathymetry option\n";
                std::abort();
            }
            
            if(bdinfo->Bot[ii].x[1] > DepthB){
                std::cout << "BELLHOP:ReadBTY: Bathymetry drops below lowest point in the sound speed profile\n";
                std::abort();
            }
        }
        
        // Convert ranges in km to m
        for(int32_t i=0; i<bdinfo->NBTYPts; ++i) bdinfo->Bot[i].x[0] *= RC(1000.0);
        
        }break;
    default:
        bdinfo->Bot = allocate<BdryPtFull>(2);
        bdinfo->Bot[0].x = vec2(-BdryInfinity(), DepthB);
        bdinfo->Bot[1].x = vec2( BdryInfinity(), DepthB);
    }
    
    ComputeBdryTangentNormal(bdinfo->Bot, false, bdinfo);
    
    if(!monotonic(&bdinfo->Bot[0].x.x, bdinfo->NBTYPts, sizeof(BdryPtFull)/sizeof(real), 0)){
        std::cout << "BELLHOP:ReadBTY: Bathymetry ranges are not monotonically increasing\n";
        std::abort();
    }
}

/**
 * Handles top and bottom boundary conditions
 * LP: Moved from readenv.cpp as it relates to boundary conditions.
 * 
 * freq: frequency [LP: :( ]
 */
inline void TopBot(const real &freq, const char (&AttenUnit)[2], real &fT, HSInfo &hs,
    LDIFile &ENVFile, std::ofstream &PRTFile, const AttenInfo *atten)
{
    real Mz, vr, alpha2_f; // values related to grain size
    
    real alphaR = RC(1500.0), betaR = RC(0.0), alphaI = RC(0.0), betaI = RC(0.0), rhoR = RC(1.0);
    real zTemp;
    
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
       std::abort();
    }
    
    // ****** Read in BC parameters depending on particular choice ******
    
    hs.cP = hs.cS = hs.rho = RC(0.0);
    
    if(hs.bc == 'A'){ // *** Half-space properties ***
        zTemp = RC(0.0);
        LIST(ENVFile); ENVFile.Read(zTemp); ENVFile.Read(alphaR);
        ENVFile.Read(betaR); ENVFile.Read(rhoR);
        ENVFile.Read(alphaI); ENVFile.Read(betaI);
        PRTFile << std::setprecision(2) << std::setw(10) << zTemp << " "
            << std::setw(10) << alphaR << " " << std::setw(10) << betaR << " "
            << std::setw(6) << rhoR << " " << std::setprecision(4) 
            << std::setw(10) << alphaI << " " << std::setw(10) << betaI << "\n";
        // dummy parameters for a layer with a general power law for attenuation
        // these are not in play because the AttenUnit for this is not allowed yet
        //freq0         = freq;
        //betaPowerLaw  = RC(1.0); //LP: Default is 1.0, this is the only other place it's set (also to 1.0).
        fT            = RC(1000.0);
        
        hs.cP  = crci(zTemp, alphaR, alphaI, freq, freq, AttenUnit, betaPowerLaw, fT, atten, PRTFile);
        hs.cS  = crci(zTemp, betaR,  betaI,  freq, freq, AttenUnit, betaPowerLaw, fT, atten, PRTFile);
        
        hs.rho = rhoR;
    }else if(hs.bc == 'G'){ // *** Grain size (formulas from UW-APL HF Handbook)
        
        // These formulas are from the UW-APL Handbook
        // The code is taken from older Matlab and is unnecesarily verbose
        // vr   is the sound speed ratio
        // rhoR is the density ratio
        LIST(ENVFile); ENVFile.Read(zTemp); ENVFile.Read(Mz);
        PRTFile << std::setprecision(2) << std::setw(10) << zTemp << " "
            << std::setw(10) << Mz << "\n";
        
        if(Mz >= RC(-1.0) && Mz < RC(1.0)){
            vr   = RC(0.002709) * SQ(Mz) - RC(0.056452) * Mz + RC(1.2778);
            rhoR = RC(0.007797) * SQ(Mz) - RC(0.17057)  * Mz + RC(2.3139);
        }else if(Mz >= RC(1.0) && Mz < RC(5.3)){
            vr   = RC(-0.0014881) * CUBE(Mz) + RC(0.0213937) * SQ(Mz) - RC(0.1382798) * Mz + RC(1.3425);
            rhoR = RC(-0.0165406) * CUBE(Mz) + RC(0.2290201) * SQ(Mz) - RC(1.1069031) * Mz + RC(3.0455);
        }else{
            vr   = RC(-0.0024324) * Mz + RC(1.0019);
            rhoR = RC(-0.0012973) * Mz + RC(1.1565);
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
        alphaR = vr * RC(1500.0);
        alphaI = alpha2_f * (vr / RC(1000.0)) * RC(1500.0) * 
            STD::log(RC(10.0)) / (RC(40.0) * M_PI); // loss parameter Sect. IV., Eq. (4) of handbook
 
        hs.cP  = crci(zTemp, alphaR, alphaI, freq, freq, {'L', ' '}, betaPowerLaw, fT, atten, PRTFile);
        hs.cS  = RC(0.0);
        hs.rho = rhoR;
    }
}
