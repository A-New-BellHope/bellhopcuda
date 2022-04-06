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
#pragma once
#include "common.hpp"
#include "attenuation.hpp"
#include "ssp.hpp"

namespace bhc {

constexpr int32_t Bdry_Number_to_Echo = 21;

template<bool THREED> struct bdry_big {};
// LP: Can't be constexpr as std::sqrt is not constexpr without GCC extensions
template<> struct bdry_big<true>  { HOST_DEVICE inline real value() { return RL(1.0e25); } };
template<> struct bdry_big<false> { HOST_DEVICE inline real value() { return STD::sqrt(REAL_MAX) / RL(1.0e5); } };
#define BDRYBIG bdry_big<THREED>::value()

/**
 * Get the top or bottom segment info (index and range interval) for range, r
 *
 * LP: t: range component of ray tangent. Endpoints of segments are handled so
 * that if the ray moves slightly along its current direction, it will remain
 * in the same segment.
 * state.lSeg: segment limits in range
 */
template<bool THREED> HOST_DEVICE inline void GetBdrySeg(real r, real t, 
    BdryStateTopBot<THREED> &bds, const BdryInfoTopBot<THREED> *bdinfotb,
    bool isTop)
{
    // LP: bdinfotb->bd.x is checked for being monotonic at load time, so we can
    // linearly search out from the last position, usually only have to move
    // by 1
    int32_t n = bdinfotb->NPts;
    bds.Iseg = bhc::max(bds.Iseg, 0);
    bds.Iseg = bhc::min(bds.Iseg, n-2);
    if(t >= FL(0.0)){
        while(bds.Iseg >= 0 && bdinfotb->bd[bds.Iseg].x.x > r) --bds.Iseg;
        while(bds.Iseg >= 0 && bds.Iseg < n-1 && bdinfotb->bd[bds.Iseg+1].x.x <= r) ++bds.Iseg;
    }else{
        while(bds.Iseg >= 0 && bds.Iseg < n-1 && bdinfotb->bd[bds.Iseg+1].x.x <  r) ++bds.Iseg;
        while(bds.Iseg >= 0 && bdinfotb->bd[bds.Iseg].x.x >= r) --bds.Iseg;
    }
    if(bds.Iseg < 0 || bds.Iseg >= n-1){
        // IsegTop MUST LIE IN [0, NATIPts-2]
        printf("Error: Get%s undefined above the ray, r=%g\n",
            isTop ? "TopSeg: Top altimetry" : "BotSeg: Bottom bathymetry", r);
        bail();
    }
    bds.lSeg.min = bdinfotb->bd[bds.Iseg  ].x.x;
    bds.lSeg.max = bdinfotb->bd[bds.Iseg+1].x.x;
    
    // LP: Only explicitly loaded in this function in 3D, loaded in containing
    // code in 2D
    bds.x = bdinfotb->bd[bds.Iseg].x.x;
    bds.n = bdinfotb->bd[bds.Iseg].x.n;
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

    Bdry[0     ].x[0] = -BDRYBIG;
    Bdry[0     ].x[1] = Bdry[1     ].x[1];
    Bdry[0     ].hs   = Bdry[1     ].hs;
    Bdry[NPts-1].x[0] =  BDRYBIG;
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
        // std::cout << "Dx, t " << Bdry[ii].Dx << " " << Bdry[ii].x << " " << (FL(1.0) / (Bdry[ii].x[1] / FL(500.0)) << "\n";
        
        // normalize the tangent vector
        Bdry[ii].Len = glm::length(Bdry[ii].t);
        Bdry[ii].t  /= Bdry[ii].Len;
        
        Bdry[ii].n[0] = (isTop ? RL(1.0) : RL(-1.0)) * Bdry[ii].t[1];
        Bdry[ii].n[1] = (isTop ? RL(-1.0) : RL(1.0)) * Bdry[ii].t[0];
    }
    
    if(CurvilinearFlag == 'C'){
        // curvilinear option: compute tangent and normal at node by averaging normals on adjacent segments
        // averaging two centered differences is equivalent to forming a single centered difference of two steps ...
        for(int32_t ii=1; ii<NPts-1; ++ii){
            sss = Bdry[ii-1].Len / (Bdry[ii-1].Len + Bdry[ii].Len);
            sss = FL(0.5); // LP: BUG? Line above is overwritten.
            Bdry[ii].Nodet = (FL(1.0) - sss) * Bdry[ii-1].t + sss * Bdry[ii].t;
        }
        
        Bdry[0     ].Nodet = vec2(FL(1.0), FL(0.0)); // tangent left-end  node
        Bdry[NPts-1].Nodet = vec2(FL(1.0), FL(0.0)); // tangent right-end node
        
        for(int32_t ii=0; ii<NPts; ++ii){
            Bdry[ii].Noden[0] = (isTop ? RL(1.0) : RL(-1.0)) * Bdry[ii].Nodet[1];
            Bdry[ii].Noden[1] = (isTop ? RL(-1.0) : RL(1.0)) * Bdry[ii].Nodet[0];
        }
        
        // compute curvature in each segment
        // LP: This allocation is not necessary, could just have two variables for
        // current and next phi. Operating on the whole array can trigger compiler
        // SIMD parallelism (AVX-512 etc.), but this is unlikely to happen for
        // atan2, and this is in one-time setup code anyway.
        phi = allocate<real>(NPts);
        // this is the angle at each node
        for(int32_t i=0; i<NPts; ++i) phi[i] = STD::atan2(Bdry[i].Nodet[1], Bdry[i].Nodet[0]);
        
        for(int32_t ii=0; ii<NPts-1; ++ii){
            Bdry[ii].kappa = (phi[ii+1] - phi[ii]) / Bdry[ii].Len; // this is curvature = dphi/ds
            Bdry[ii].Dxx   = (Bdry[ii+1].Dx - Bdry[ii].Dx) / // second derivative
                             (Bdry[ii+1].x[0] - Bdry[ii].x[0]); 
            Bdry[ii].Dss   = Bdry[ii].Dxx * CUBE(Bdry[ii].t[0]); // derivative in direction of tangent
            //std::cout << "kappa, Dss, Dxx " << Bdry[ii].kappa << " " << Bdry[ii].Dss << " " << Bdry[ii].Dxx
            //    << " " << FL(1.0) / ((FL(8.0) / SQ(FL(1000.0))) * CUBE(STD::abs(Bdry[ii].x[1])))
            //    << " " << Bdry[ii].x[1] << " "
            //    << FL(-1.0) / (FL(4.0) * CUBE(Bdry[ii].x[1]) / FL(1000000.0))
            //    << " " << Bdry[ii].x[1] << "\n";
            
            Bdry[ii].kappa = Bdry[ii].Dss; // over-ride kappa !!!!!
        }
        
        deallocate(phi);
    }else{
        for(int32_t i=0; i<NPts; ++i) Bdry[i].kappa = FL(0.0);
    }
}

template<bool THREED> inline void ReadBoundary(std::string FileRoot, char BdryDefMode, real BdryDepth,
    PrintFileEmu &PRTFile, BdryInfoTopBot<THREED> *bdinfotb, bool isTop)
{
    const char *s_atibty = isTop ? "ati" : "bty";
    const char *s_ATIBTY = isTop ? "ATI" : "BTY";
    const char *s_altimetrybathymetry = isTop ? "altimetry" : "bathymetry";
    const char *s_AltimetryBathymetry = isTop ? "Altimetry" : "Bathymetry";
    const char *s_topbottom = isTop ? "top" : "bottom";
    
    switch(BdryDefMode){
    case '~':
    case '*':{
        if constexpr(THREED){
            PRTFile << "*********************************\n";
        }else{
            PRTFile << "__________________________________________________________________________\n\n";
        }
        PRTFile << "Using " << s_topbottom << "-" << s_altimetrybathymetry << " file\n";
        
        LDIFile BDRYFile(FileRoot + "." + s_atibty);
        if(!BDRYFile.Good()){
            PRTFile << s_ATIBTY << "File = " << FileRoot << "." << s_atibty << "\n";
            std::cout << "Read" << s_ATIBTY << ": Unable to open " 
                << s_altimetrybathymetry << " file\n";
            std::abort();
        }
        
        LIST(BDRYFile); BDRYFile.Read(bdinfotb->type, THREED ? 1 : 2);
        if constexpr(THREED) bdinfotb->type[1] = ' ';
        switch(bdinfotb->type[0]){
        case 'R':
            if constexpr(THREED){
                PRTFile << "Regular grid for a 3D run\n";
            }else{
                PRTFile << s_atibty << "Type R not supported for 2D runs\n"; std::abort();
            }
            break;
        case 'C':
            if constexpr(THREED){
                PRTFile << "Regular grid for a 3D run (curvilinear)\n";
            }else{
                PRTFile << "Curvilinear Interpolation\n";
            }
            break;
        case 'L':
            if constexpr(THREED){
                PRTFile << s_atibty << "Type L not supported for 3D runs\n"; std::abort();
            }else{
                PRTFile << "Piecewise linear interpolation\n";
            }
            break;
        default:
            std::cout << "Read" << s_ATIBTY << ": Unknown option for selecting " 
                << s_altimetrybathymetry << " interpolation\n";
            std::abort();
        }
        
        if constexpr(THREED){
            // x values
            LIST(BDRYFile); BDRYFile.Read(bdinfotb->NPts.x);
            PRTFile << "\nNumber of " << s_altimetrybathymetry << " points in x-direction "
                << bdinfotb->NPts.x << "\n";
            
            real *Globalx = allocate<real>(std::max(bdinfotb->NPts.x, 3));
            Globalx[2] = FL(-999.9);
            LIST(BDRYFile); BDRYFile.Read(Globalx, bdinfotb->NPts.x);
            SubTab(Globalx, bdinfotb->NPts.x);
            EchoVector(Globalx, bdinfotb->NPts.x, PRTFile, Bdry_Number_to_Echo);
            
            // y values
            LIST(BDRYFile); BDRYFile.Read(bdinfotb->NPts.y);
            PRTFile << "\nNumber of " << s_altimetrybathymetry << " points in y-direction "
                << bdinfotb->NPts.y << "\n";
            
            real *Globaly = allocate<real>(std::max(bdinfotb->NPts.y, 3));
            Globaly[2] = FL(-999.9);
            LIST(BDRYFile); BDRYFile.Read(Globaly, bdinfotb->NPts.y);
            SubTab(Globalx, bdinfotb->NPts.y);
            EchoVector(Globalx, bdinfotb->NPts.y, PRTFile, Bdry_Number_to_Echo);
            
            // convert km to m
            for(int32_t i=0; i<bdinfotb->NPts.x; ++i) Globalx[i] *= FL(1000.0);
            for(int32_t i=0; i<bdinfotb->NPts.y; ++i) Globaly[i] *= FL(1000.0);
            
            // z values
            checkallocate(bdinfotb->bd, bdinfotb->NPts.x * bdinfotb->NPts.y);
            
            PRTFile << "\n";
            bool warnedNaN = false;
            for(int32_t iy=0; iy<bdinfotb->NPts.y; ++iy){
                LIST(BDRYFile); // read a row of depths
                for(int32_t ix=0; ix<bdinfotb->NPts.x; ++ix){
                    vec3 &x = bdinfotb->bd[ix*bdinfotb->NPts.y+iy].x;
                    BDRYFile.Read(x.z);
                    if(!std::isfinite(x.z) && !warnedNaN){
                        PRTFile << "Warning in " BHC_PROGRAMNAME "3D - Read" << s_ATIBTY 
                            << "3D : The " << s_altimetrybathymetry << " file contains a NaN\n";
                        warnedNaN = true;
                    }
                    x.x = Globalx[ix];
                    x.y = Globaly[iy];
                }
            }
        }else{
        
            LIST(BDRYFile); BDRYFile.Read(bdinfotb->NPts);
            PRTFile << "Number of " << s_altimetrybathymetry << " points = " 
                << bdinfotb->NPts << "\n";
            bdinfotb->NPts += 2; // we'll be extending the s_altimetrybathymetry to infinity to the left and right
            
            checkallocate(bdinfotb->bd, bdinfotb->NPts);
            
            // LP: BUG: Geoacoustics are supported for altimetry, but the
            // header for geoacoustics is only supported for bathymetry.
            if(isTop || bdinfotb->type[1] == 'S' || bdinfotb->type[1] == ' '){
                if(!isTop){
                    PRTFILE << "Short format (" << s_altimetrybathymetry << " only)\n";
                }
                PRTFile << "\n Range (km)  Depth (m)\n";
            }else if(bdinfotb->type[1] == 'L'){
                PRTFile << "Long format (" << s_altimetrybathymetry << " and geoacoustics)\n";
                PRTFile << "Range (km)  Depth (m)  alphaR (m/s)  betaR  rho (g/cm^3)  alphaI     betaI\n";
            }else{
                std::cout << "Read" << s_ATIBTY << ": Unknown option for selecting " 
                    << s_altimetrybathymetry << " option\n";
                std::abort();
            }
            
            for(int32_t ii=1; ii<bdinfotb->NPts-1; ++ii){
                switch(bdinfotb->type[1]){
                case 'S':
                case ' ':
                    LIST(BDRYFile); BDRYFile.Read(bdinfotb->bd[ii].x);
                    // LP: This condition was previously ii == bdinfotb->NPts - 1,
                    // which will never be satisfied due to the loop bounds
                    if(ii < Bdry_Number_to_Echo || ii == bdinfotb->NPts - 2){ // echo some values
                        PRTFile << std::setprecision(3) << bdinfotb->bd[ii].x << "\n";
                    }
                    break;
                case 'L':
                    LIST(BDRYFile); BDRYFile.Read(bdinfotb->bd[ii].x);
                    BDRYFile.Read(bdinfotb->bd[ii].hs.alphaR);
                    BDRYFile.Read(bdinfotb->bd[ii].hs.betaR);
                    BDRYFile.Read(bdinfotb->bd[ii].hs.rho);
                    BDRYFile.Read(bdinfotb->bd[ii].hs.alphaI);
                    BDRYFile.Read(bdinfotb->bd[ii].hs.betaI);
                    // LP: Same change as above
                    if(ii < Bdry_Number_to_Echo || ii == bdinfotb->NPts - 2){ // echo some values
                        PRTFile << std::setprecision(3) << bdinfotb->bd[ii].x << " "
                            << bdinfotb->bd[ii].hs.alphaR << " "
                            << bdinfotb->bd[ii].hs.betaR << " "
                            << bdinfotb->bd[ii].hs.rho << " "
                            << bdinfotb->bd[ii].hs.alphaI << " "
                            << bdinfotb->bd[ii].hs.betaI << "\n";
                    }
                    break;
                default:
                    std::cout << "Read" << s_ATIBTY << ": Unknown option for selecting " 
                        << s_altimetrybathymetry << " option\n";
                    std::abort();
                }
                
                if(bdinfotb->bd[ii].x[1] < BdryDepth){
                    std::cout << "BELLHOP:Read" << s_ATIBTY << ": " << s_AltimetryBathymetry 
                        << " rises above highest point in the sound speed profile\n";
                    std::abort();
                }
            }
            
            if(!monotonic(&bdinfotb->bd[0].x.x, bdinfotb->NPts, sizeof(BdryPtFull)/sizeof(real), 0)){
                std::cout << "BELLHOP:Read" << s_ATIBTY << ": " << s_AltimetryBathymetry 
                    << " ranges are not monotonically increasing\n";
                std::abort();
            }
            
            // Convert ranges in km to m
            for(int32_t i=0; i<bdinfotb->NPts; ++i) bdinfotb->bd[i].x[0] *= FL(1000.0);
            
        }
        
        }break;
    default:
        if constexpr(THREED){
            bdinfotb->atiType[0] = 'R';
            bdinfotb->NPts = int2(2, 2);
            checkallocate(bdinfotb->bd, 2*2);
            
            // LP: TODO/BUG: Top_deltax and Top_deltay initialized here. This
            // value is only used if the ray goes outside the region where
            // altimetry is defined. But there's 2 problems with this:
            // 1) It will contain the leftover value from a previous step--or a
            //    previous ray!--rather than the initial value from here.
            // 2) This initial value is only written on this codepath with no
            //    altimetry data--it's not written on the other codepath where
            //    altimetry is defined. This means, on that codepath, if the
            //    initial value would ever be used, it will be uninitialized.
            //    Furthermore, it is much more likely the ray will go out of
            //    bounds of the altimetry definition when it is manually defined
            //    over a specific region, than on this codepath where the
            //    altimetry is defined to be 10,000x the diameter of the Milky
            //    Way!
            //Top_deltax = FL(2.0) * BDRYBIG;
            //Top_deltay = FL(2.0) * BDRYBIG;
            
            bdinfotb->bd[0].x = vec3(-BDRYBIG, -BDRYBIG, BdryDepth);
            bdinfotb->bd[1].x = vec3(-BDRYBIG,  BDRYBIG, BdryDepth);
            bdinfotb->bd[2].x = vec3( BDRYBIG, -BDRYBIG, BdryDepth);
            bdinfotb->bd[3].x = vec3( BDRYBIG,  BDRYBIG, BdryDepth);
            
            for(int32_t i=0; i<4; ++i){
                bdinfotb->bd[i].t  = vec3(FL(1.0), FL(0.0), FL( 0.0));
                bdinfotb->bd[i].n1 = vec3(FL(0.0), FL(0.0), FL(-1.0));
                bdinfotb->bd[i].n2 = vec3(FL(0.0), FL(0.0), FL(-1.0));
            }
            
            return; // LP: No ComputeBdryTangentNormal cause done manually here
        }else{
            checkallocate(bdinfotb->bd, 2);
            bdinfotb->bd[0].x = vec2(-BDRYBIG, BdryDepth);
            bdinfotb->bd[1].x = vec2( BDRYBIG, BdryDepth);
        }
    }
    
    ComputeBdryTangentNormal(bdinfotb->bd, isTop, bdinfotb);
    
    // LP: TODO/BUG: 3D version has initialization for xTopSeg / yTopSeg here,
    // which probably also means state is carried over from one ray to the next
}

/**
 * Handles top and bottom boundary conditions
 * LP: Moved from readenv.cpp as it relates to boundary conditions.
 * 
 * freq: center / nominal frequency (wideband not supported)
 */
inline void TopBot(const real &freq, const char (&AttenUnit)[2], real &fT, HSInfo &hs,
    LDIFile &ENVFile, PrintFileEmu &PRTFile, const AttenInfo *atten,  HSInfo &RecycledHS)
{
    real Mz, vr, alpha2_f; // values related to grain size
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
        PRTFile << "    Writing an IFL file\n"; break;
    case 'P':
        PRTFile << "    reading PRECALCULATED IFL\n"; break;
    default:
       std::cout << "TopBot: Unknown boundary condition type\n";
       std::abort();
    }
    
    // ****** Read in BC parameters depending on particular choice ******
    
    hs.cP = hs.cS = hs.rho = FL(0.0);
    
    if(hs.bc == 'A'){ // *** Half-space properties ***
        zTemp = FL(0.0);
        LIST(ENVFile); ENVFile.Read(zTemp); ENVFile.Read(RecycledHS.alphaR);
        ENVFile.Read(RecycledHS.betaR); ENVFile.Read(RecycledHS.rho);
        ENVFile.Read(RecycledHS.alphaI); ENVFile.Read(RecycledHS.betaI);
        PRTFile << std::setprecision(2) << std::setw(10) << zTemp << " "
            << std::setw(10) << RecycledHS.alphaR << " " << std::setw(10) << RecycledHS.betaR << " "
            << std::setw(6) << RecycledHS.rho << " " << std::setprecision(4) 
            << std::setw(10) << RecycledHS.alphaI << " " << std::setw(10) << RecycledHS.betaI << "\n";
        // dummy parameters for a layer with a general power law for attenuation
        // these are not in play because the AttenUnit for this is not allowed yet
        //freq0         = freq;
        //betaPowerLaw  = FL(1.0); //LP: Default is 1.0, this is the only other place it's set (also to 1.0).
        fT            = FL(1000.0);
        
        hs.cP  = crci(zTemp, RecycledHS.alphaR, RecycledHS.alphaI, freq, freq, AttenUnit, betaPowerLaw, fT, atten, PRTFile);
        hs.cS  = crci(zTemp, RecycledHS.betaR,  RecycledHS.betaI,  freq, freq, AttenUnit, betaPowerLaw, fT, atten, PRTFile);
        // printf("%g %g %g %g %c%c %g %g\n", zTemp, RecycledHS.alphaR, RecycledHS.alphaI, freq,
        //     AttenUnit[0], AttenUnit[1], betaPowerLaw, fT);
        // printf("cp computed to (%g,%g)\n", hs.cP.real(), hs.cP.imag());
        
        hs.rho = RecycledHS.rho;
    }else if(hs.bc == 'G'){ // *** Grain size (formulas from UW-APL HF Handbook)
        
        // These formulas are from the UW-APL Handbook
        // The code is taken from older Matlab and is unnecesarily verbose
        // vr   is the sound speed ratio
        // rho is the density ratio
        LIST(ENVFile); ENVFile.Read(zTemp); ENVFile.Read(Mz);
        PRTFile << std::setprecision(2) << std::setw(10) << zTemp << " "
            << std::setw(10) << Mz << "\n";
        
        if(Mz >= FL(-1.0) && Mz < FL(1.0)){
            vr             = FL(0.002709) * SQ(Mz) - FL(0.056452) * Mz + FL(1.2778);
            RecycledHS.rho = FL(0.007797) * SQ(Mz) - FL(0.17057)  * Mz + FL(2.3139);
        }else if(Mz >= FL(1.0) && Mz < FL(5.3)){
            vr             = FL(-0.0014881) * CUBE(Mz) + FL(0.0213937) * SQ(Mz) - FL(0.1382798) * Mz + FL(1.3425);
            RecycledHS.rho = FL(-0.0165406) * CUBE(Mz) + FL(0.2290201) * SQ(Mz) - FL(1.1069031) * Mz + FL(3.0455);
        }else{
            vr             = FL(-0.0024324) * Mz + FL(1.0019);
            RecycledHS.rho = FL(-0.0012973) * Mz + FL(1.1565);
        }
        
        if(Mz >= FL(-1.0) && Mz < FL(0.0)){
            alpha2_f = FL(0.4556);
        }else if(Mz >= FL(0.0) && Mz < FL(2.6)){
            alpha2_f = FL(0.4556) + FL(0.0245) * Mz;
        }else if(Mz >= FL(2.6) && Mz < FL(4.5)){
            alpha2_f = FL(0.1978) + FL(0.1245) * Mz;
        }else if(Mz >= FL(4.5) && Mz < FL(6.0)){
            alpha2_f = FL(8.0399) - FL(2.5228) * Mz + FL(0.20098) * SQ(Mz);
        }else if(Mz >= FL(6.0) && Mz < FL(9.5)){
            alpha2_f = FL(0.9431) - FL(0.2041) * Mz + FL(0.0117) * SQ(Mz);
        }else{
            alpha2_f =  FL(0.0601);
        }
        
        // AttenUnit = 'L';  // loss parameter
        // !! following uses a reference sound speed of 1500 ???
        // !! should be sound speed in the water, just above the sediment
        // the term vr / 1000 converts vr to units of m per ms 
        RecycledHS.alphaR = vr * FL(1500.0);
        RecycledHS.alphaI = alpha2_f * (vr / FL(1000.0)) * FL(1500.0) * 
            STD::log(FL(10.0)) / (FL(40.0) * REAL_PI); // loss parameter Sect. IV., Eq. (4) of handbook
 
        hs.cP  = crci(zTemp, RecycledHS.alphaR, RecycledHS.alphaI, freq, freq, {'L', ' '}, betaPowerLaw, fT, atten, PRTFile);
        hs.cS  = FL(0.0);
        hs.rho = RecycledHS.rho;
    }
}

}
