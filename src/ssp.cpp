#include "ssp.hpp"

#define READ_SSP_ARGS real Depth, real freq, const real &fT, SSPStructure *ssp, \
    LDIFile &ENVFile, std::ofstream &PRTFile, const AttenInfo *atten
#define CALL_READ_SSP_ARGS Depth, freqinfo->freq0, fT, ssp, ENVFile, PRTFile, atten

/**
 * reads the SSP data from the environmental file and convert to Nepers/m
 */
void ReadSSP(READ_SSP_ARGS)
{
    PRTFile << "\nSound speed profile:\n";
    PRTFile << "   z (m)     alphaR (m/s)   betaR  rho (g/cm^3)  alphaI     betaI\n";
    
    ssp->NPts = 1;
    
    for(int32_t iz=0; iz<MaxSSP; ++iz){
        real alphaR, betaR, rhoR, alphaI, betaI;
        ENVFile.List(); ENVFile.Read(ssp->z[iz]); 
        ENVFile.Read(alphaR); ENVFile.Read(betaR); ENVFile.Read(rhoR);
        ENVFile.Read(alphaI); ENVFile.Read(betaI);
        PRTFile << std::setprecision(2) << ssp->z[iz] << " " << alphaR << " " << betaR
            << " " << rhoR << " " << std::setprecision(4) << alphaI << " " << betaI << "\n";
        
        ssp->c[iz] = crci(ssp->z[iz], alphaR, alphaI, freq, freq, ssp->AttenUnit,
            betaPowerLaw, fT, atten, PRTFile);
        ssp->rho[iz] = rhoR;
        
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
        if(std::abs(ssp->z[iz] - Depth) < RC(100.0) * REAL_EPSILON){
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
    real Depth = x[1];
    ReadSSP(CALL_READ_SSP_ARGS);
    
    for(int32_t i=0; i<ssp->NPts; ++i) ssp->n2[i] = RC(1.0) / SQ(ssp->c[i]);
    
    // compute gradient, n2z
    for(int32_t iz=1; iz<ssp->NPts; ++iz){
        ssp->n2z[iz-1] = (ssp->n2[iz] - ssp->n2[iz-1]) /
                         (ssp->z [iz] - ssp->z [iz-1]);
    }
}

void InitcLinear(SSP_INIT_ARGS)
{
    real Depth = x[1];
    ReadSSP(CALL_READ_SSP_ARGS);
}

void InitcPCHIP(SSP_INIT_ARGS)
{
    real Depth = x[1];
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
    real Depth = x[1];
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
    real Depth = x[1];
    ReadSSP(CALL_READ_SSP_ARGS);
    
    // Read the 2D SSP matrix
    PRTFile << "__________________________________________________________________________\n\n";
    PRTFile << "Using range-dependent sound speed\n";
    
    LDIFile SSPFile(FileRoot + ".ssp");
    SSPFile.List(); SSPFile.Read(ssp->Nr);
    PRTFile << "Number of SSP ranges = " << ssp->Nr << "\n";
    
    if(ssp->Nr < 2){
        PRTFile << "READIN: Quad: You must have a least two profiles in your 2D SSP field\n";
        std::abort();
    }
    
    ssp->cMat  = allocate<real>( ssp->NPts    * ssp->Nr);
    ssp->czMat = allocate<real>((ssp->NPts-1) * ssp->Nr);
    ssp->Seg.r = allocate<real>(ssp->Nr);
    
    SSPFile.List(); SSPFile.Read(ssp->Seg.r, ssp->Nr);
    PRTFile << "\nProfile ranges (km):\n" << std::setprecision(2);
    for(int32_t i=0; i<ssp->Nr; ++i) PRTFile << ssp->Seg.r[i] << " ";
    PRTFile << "\n";
    
    for(int32_t i=0; i<ssp->Nr; ++i) ssp->Seg.r[i] *= RC(1000.0); // convert km to m
    
    PRTFile << "\nSound speed matrix:\n";
    PRTFile << " Depth (m )     Soundspeed (m/s)\n";
    for(int32_t iz2=0; iz2<ssp->NPts; ++iz2){
        SSPFile.List(); SSPFile.Read(&ssp->cMat[iz2*ssp->Nr], ssp->Nr);
        // PRTFile << "iSegz depth = " << std::setprecision(2) << ssp->z[iz2] << " m\n";
        PRTFile << std::setprecision(2) << ssp->z[iz2] << " ";
        for(int32_t i=0; i<ssp->NPts; ++i) PRTFile << ssp->cMat[iz2*ssp->Nr+i] << " ";
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

void InitializeSSP(SSP_INIT_ARGS)
{
    switch(ssp->Type){
    case 'N': // N2-linear profile option
        Initn2Linear(SSP_CALL_INIT_ARGS); break;
    case 'C': // C-linear profile option
        InitcLinear (SSP_CALL_INIT_ARGS); break;
    case 'P': // monotone PCHIP ACS profile option
        InitcPCHIP  (SSP_CALL_INIT_ARGS); break;
    case 'S': // Cubic spline profile option
        InitcCubic  (SSP_CALL_INIT_ARGS); break;
    case 'Q':
        InitQuad    (SSP_CALL_INIT_ARGS); break;
    /* case 'H':
        // this is called by BELLHOP3D only once, during READIN
        // possibly the logic should be changed to call EvaluateSSP2D or 3D
        x3 = vec3(RC(0.0), RC(0.0), x.y);
        InitHexahedral(x3, freq, ssp); break; */
    case 'A': // Analytic profile option
        break; //LP: No init for analytic.
    default:
        printf("InitializeSSP: Invalid profile option %c\n", ssp->Type);
        std::abort();
    }
}
