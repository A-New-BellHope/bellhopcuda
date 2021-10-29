#pragma once
#include "common.hpp"
#include "ldio.hpp"
#include "ssp.hpp"

struct BdryPtFull {
    vec2 x, t, n; // coordinate, tangent, and outward normal for a segment
    vec2 Nodet, Noden; // tangent and normal at the node, if the curvilinear option is used
    real Len, Kappa; // length and curvature of a segement
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
    const BdryInfo *bdry)
{
    // LP: bdry->Top.x is checked for being monotonic at load time, so we can
    // linearly search out from the last position, usually only have to move
    // by 1
    int32_t n = bdry->NATIPts;
    IsegTop = STD::max(IsegTop, 0);
    IsegTop = STD::min(IsegTop, n-2);
    while(IsegTop >= 0 && bdry->Top[IsegTop].x.x > r) --IsegTop;
    while(IsegTop >= 0 && IsegTop < n-1 && bdry->Top[IsegTop+1].x.x < r) ++IsegTop;
    if(IsegTop < 0 || IsegTop >= n-1){
        // IsegTop MUST LIE IN [0, NatiPts-2]
        printf("Error: GetTopSeg: Top altimetry undefined above the ray, r=%f\n", r);
        bail();
    }
    rTopSeg.x = bdry->Top[IsegTop].x.x;
    rTopSeg.y = bdry->Top[IsegTop+1].x.x;
}

/**
 * Get the Bottom segment info (index and range interval) for range, r
 *
 * rBotSeg: segment limits in range
 */
HOST_DEVICE inline void GetBotSeg(real r, int32_t &IsegBot, vec2 &rBotSeg,
    const BdryInfo *bdry)
{
    // LP: bdry->Bot.x is checked for being monotonic at load time, so we can
    // linearly search out from the last position, usually only have to move
    // by 1
    int32_t n = bdry->NATIPts;
    IsegBot = STD::max(IsegBot, 0);
    IsegBot = STD::min(IsegBot, n-2);
    while(IsegBot >= 0 && bdry->Bot[IsegBot].x.x > r) --IsegBot;
    while(IsegBot >= 0 && IsegBot < n-1 && bdry->Bot[IsegBot+1].x.x < r) ++IsegBot;
    if(IsegBot < 0 || IsegBot >= n-1){
        // IsegBot MUST LIE IN [0, NatiPts-2]
        printf("Error: GetBotSeg: Bottom altimetry undefined below the source, r=%f\n", r);
        bail();
    }
    rBotSeg.x = bdry->Bot[IsegBot].x.x;
    rBotSeg.y = bdry->Bot[IsegBot+1].x.x;
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
    
    for(int32_t ii=0; ii<Npts-1; ++ii){
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
        
        Bdry[ii].Noden[0] = (isTop ? RC(1.0) : RC(-1.0)) * Bdry[ii].Nodet[1];
        Bdry[ii].Noden[1] = (isTop ? RC(-1.0) : RC(1.0)) * Bdry[ii].Nodet[0];
        
        // compute curvature in each segment
        phi = allocate<real>(Npts);
        // this is the angle at each node
        for(int32_t i=0; i<Npts; ++i) phi[i] = STD::atan2(Bdry[i].Nodet[1], Bdry[i].Nodet[0]);
        
        for(int32_t ii=0; ii<Npts-1; ++ii){
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
        for(int32_t i=0; i<Npts; ++i) Bdry[i].kappa = RC(0.0);
    }
}

inline void ReadATI(std::string FileRoot, char TopATI, real DepthT,
    std::ofstream &PRTFile, BdryInfo *bdry)
{
    //LP: Removed phi, which got allocated but not used. Just to check that
    //there was enough memory for it to be allocated later? Or bug?
    switch(TopATI){
    case '~':
    case '*':
        PRTFile << "__________________________________________________________________________\n\n";
        PRTFile << "Using top-altimetry file\n";
        
        LDIFile ATIFile(FileRoot + ".ati");
        if(!ATIFile.Good()){
            PRTFile << "ATIFile = " << FileRoot << ".ati\n";
            std::cout << "ReadATI: Unable to open altimetry file\n";
            std::abort();
        }
        
        ATIFile.List(); ATIFile.Read(bdry->atiType, 2);
        switch(bdry->atiType[0]){
        case 'C':
            PRTFile << "Curvilinear Interpolation\n"; break;
        case 'L':
            PRTFile << "Piecewise linear interpolation\n"; break;
        default:
            std::cout << "ReadATI: Unknown option for selecting altimetry interpolation\n";
            std::abort();
        }
        
        ATIFile.List(); ATIFile.Read(bdry->NATIPts);
        PRTFile << "Number of altimetry points = " << bdry->NATIPts << "\n";
        bdry->NATIPts += 2; // we'll be extending the altimetry to infinity to the left and right
        
        bdry->Top = allocate<BdryPtFull>(bdry->NATIPts);
        
        PRTFile << "\n Range (km)  Depth (m)\n";
        
        for(int32_t ii=1; ii<bdry->NATIPts-1; ++ii){
            switch(atiType[1]){
            case 'S':
            case '\0':
                ATIFile.List(); ATIFile.Read(bdry->Top[ii].x);
                if(ii < Bdry_Number_to_Echo || ii == bdry->NATIPts){ // echo some values
                    //LP: Bug? Even in the Fortran, ii should never equal
                    //NATIPts as the loop ends at NATIPts-1.
                    PRTFile << std::setprecision(3) << bdry->Top[ii].x << "\n";
                }
                break;
            case 'L':
                ATIFile.List(); ATIFile.Read(bdry->Top[ii].x);
                ATIFile.Read(bdry->Top[ii].hs.alphaR);
                ATIFile.Read(bdry->Top[ii].hs.betaR);
                ATIFile.Read(bdry->Top[ii].hs.rho);
                ATIFile.Read(bdry->Top[ii].hs.alphaI);
                ATIFile.Read(bdry->Top[ii].hs.betaI);
                if(ii < Bdry_Number_to_Echo || ii == bdry->NATIPts){ // echo some values
                    PRTFile << std::setprecision(3) << bdry->Top[ii].x << " "
                        << bdry->Top[ii].hs.alphaR << " "
                        << bdry->Top[ii].hs.betaR << " "
                        << bdry->Top[ii].hs.rho << " "
                        << bdry->Top[ii].hs.alphaI << " "
                        << bdry->Top[ii].hs.betaI << "\n";
                }
                break;
            default:
                std::cout << "ReadATI: Unknown option for selecting altimetry option\n";
                std::abort();
            }
            
            if(bdry->Top[ii].x[1] < DepthT){
                std::cout << "BELLHOP:ReadATI: Altimetry rises above highest point in the sound speed profile\n";
                std::abort();
            }
        }
        
        // Convert ranges in km to m
        for(int32_t i=0; i<bdry->NATIPts; ++i) bdry->Top[i].x[0] *= RC(1000.0);
        
        break;
    default:
        bdry->Top = allocate<BdryPtFull>(2);
        bdry->Top[0].x = vec2(-STD::sqrt(REAL_MAX) / RC(1.0e5), DepthT);
        bdry->Top[0].x = vec2( STD::sqrt(REAL_MAX) / RC(1.0e5), DepthT);
    }
    
    ComputeBdryTangentNormal(bdry->Top, true);
    
    if(!monotonic(bdry->Top.x, bdry->NATIPts, 0)){
        std::cout << "BELLHOP:ReadATI: Altimetry ranges are not monotonically increasing\n";
        std::abort();
    }
}

inline void ReadBTY(std::string FileRoot, char BotBTY, real DepthB,
    std::ofstream &PRTFile, BdryInfo *bdry)
{
    switch(BotBTY){
    case '~':
    case '*':
        PRTFile << "__________________________________________________________________________\n\n";
        PRTFile << "Using bottom-bathymetry file\n";
        
        LDIFile BTYFile(FileRoot + ".bty");
        if(!BTYFile.Good()){
            PRTFile << "BTYFile = " << FileRoot << ".bty\n";
            std::cout << "ReadATI: Unable to open bathymetry file\n";
            std::abort();
        }
        
        BTYFile.List(); BTYFile.Read(bdry->btyType, 2);
        
        switch(bdry->btyType[0]){
        case 'C':
            PRTFile << "Curvilinear Interpolation\n"; break;
        case 'L':
            PRTFile << "Piecewise linear interpolation\n"; break;
        default:
            std::cout << "ReadBTY: Unknown option for selecting bathymetry interpolation\n";
            std::abort();
        }
        
        BTYFile.List(); BTYFile.Read(bdry->NBTYPts);
        PRTFile << "Number of bathymetry points = " << bdry->NBTYPts << "\n";
        
        bdry->NBTYPts += 2; // we'll be extending the bathymetry to infinity to the left and right
        bdry->Bot = allocate<BdryPtFull>(bdry->NBTYPts);
        
        PRTFile << "\n";
        switch(btyType[1]){
        case 'S':
        case '\0':
            PRTFile << "Short format (bathymetry only)\n";
            PRTFile << " Range (km)  Depth (m)\n";
            break;
        case 'L':
            PRTFile << "Long format (bathymetry and geoacoustics)\n";
            PRTFile << "Range (km)  Depth (m)  alphaR (m/s)  betaR  rho (g/cm^3)  alphaI     betaI\n";
            break;
        default:
            std::cout << "ReadBTY: Unknown option for selecting bathymetry interpolation";
            std::abort();
        }
        
        for(int32_t ii=1; ii<bdry->NBTYPts-1; ++ii){
            switch(atiType[1]){
            case 'S':
            case '\0':
                BTYFile.List(); BTYFile.Read(bdry->Bot[ii].x);
                if(ii < Bdry_Number_to_Echo || ii == bdry->NBTYPts){ // echo some values
                    //LP: Bug? Even in the Fortran, ii should never equal
                    //NBTYPts as the loop ends at NBTYPts-1.
                    PRTFile << std::setprecision(3) << bdry->Bot[ii].x << "\n";
                }
                break;
            case 'L':
                BTYFile.List(); BTYFile.Read(bdry->Bot[ii].x);
                BTYFile.Read(bdry->Bot[ii].hs.alphaR);
                BTYFile.Read(bdry->Bot[ii].hs.betaR);
                BTYFile.Read(bdry->Bot[ii].hs.rho);
                BTYFile.Read(bdry->Bot[ii].hs.alphaI);
                BTYFile.Read(bdry->Bot[ii].hs.betaI);
                if(ii < Bdry_Number_to_Echo || ii == bdry->NBTYPts){ // echo some values
                    PRTFile << std::setprecision(2) << bdry->Bot[ii].x << " "
                        << bdry->Bot[ii].hs.alphaR << " "
                        << bdry->Bot[ii].hs.betaR << " "
                        << bdry->Bot[ii].hs.rho << " " << std::setprecision(4)
                        << bdry->Bot[ii].hs.alphaI << " "
                        << bdry->Bot[ii].hs.betaI << "\n";
                }
                break;
            default:
                std::cout << "ReadBTY: Unknown option for selecting bathymetry option\n";
                std::abort();
            }
            
            if(bdry->Bot[ii].x[1] > DepthB){
                std::cout << "BELLHOP:ReadBTY: Bathymetry drops below lowest point in the sound speed profile\n";
                std::abort();
            }
        }
        
        // Convert ranges in km to m
        for(int32_t i=0; i<bdry->NBTYPts; ++i) bdry->Bot[i].x[0] *= RC(1000.0);
        
        break;
    default:
        bdry->Bot = allocate<BdryPtFull>(2);
        bdry->Bot[0].x = vec2(-BdryInfinity(), DepthB);
        bdry->Bot[0].x = vec2( BdryInfinity(), DepthB);
    }
    
    ComputeBdryTangentNormal(bdry->Bot, false);
    
    if(!monotonic(bdry->Bot.x, bdry->NBTYPts, 0)){
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
    LDIFile &ENVFile, std::ofstream &PRTFile)
{
    real Mz, vr, alpha2_f; // values related to grain size
    
    real alphaR = RC(1500.0), betaR = RC(0.0), alphaI = RC(0.0), betaI = RC(0.0), rhoR = RC(1.0);
    real ztemp;
    
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
        zTemp = RC(0.0);
        ENVFile.List(); ENVFile.Read(zTemp); ENVFile.Read(alphaR);
        ENVFile.Read(betaR); ENVFile.read(rhoR);
        ENVFile.read(alphaI); ENVFile.read(beta);
        PRTFile << std::setprecision(2) << std::setw(10) << ztemp << " "
            << std::setw(10) << alphaR << " " << std::setw(10) << betaR << " "
            << std::setw(6) << rhoR << " " << std::setprecision(4) 
            << std::setw(10) << alphaI << " " << std::setw(10) << betaI << "\n";
        // dummy parameters for a layer with a general power law for attenuation
        // these are not in play because the AttenUnit for this is not allowed yet
        //freq0         = freq;
        //betaPowerLaw  = RC(1.0); //LP: Default is 1.0, this is the only other place it's set (also to 1.0).
        fT            = RC(1000.0);
        
        hs.cp  = crci(zTemp, alphaR, alphaI, freq, freq, AttenUnit, betaPowerLaw, fT);
        hs.cs  = crci(zTemp, betaR,  betaI,  freq, freq, AttenUnit, betaPowerLaw, fT);
        
        hs.rho = rhoR;
    }else if(hs.bc == 'G'){ // *** Grain size (formulas from UW-APL HF Handbook)
        
        // These formulas are from the UW-APL Handbook
        // The code is taken from older Matlab and is unnecesarily verbose
        // vr   is the sound speed ratio
        // rhor is the density ratio
        ENVFile.List(); ENVFile.Read(zTemp); ENVFile.Read(Mz);
        PRTFile << std::setprecision(2) << std::setw(10) << zTemp << " "
            << std::setw(10) << Mz << "\n";
        
        if(Mz >= RC(-1.0) && Mz < RC(1.0)){
            vr   = RC(0.002709) * SQ(Mz) - RC(0.056452) * Mz + RC(1.2778);
            rhor = RC(0.007797) * SQ(Mz) - RC(0.17057)  * Mz + RC(2.3139);
        }else if(Mz >= RC(1.0) && Mz < RC(5.3)){
            vr   = RC(-0.0014881) * Cube(Mz) + RC(0.0213937) * SQ(Mz) - RC(0.1382798) * Mz + RC(1.3425);
            rhor = RC(-0.0165406) * Cube(Mz) + RC(0.2290201) * SQ(Mz) - RC(1.1069031) * Mz + RC(3.0455);
        }else{
            vr   = RC(-0.0024324) * Mz + RC(1.0019);
            rhor = RC(-0.0012973) * Mz + RC(1.1565);
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
 
        hs.cp  = crci(zTemp, alphaR, alphaI, freq, freq, "L ", betaPowerLaw, fT);
        hs.cs  = RC(0.0);
        hs.rho = rhoR;
    }
}
