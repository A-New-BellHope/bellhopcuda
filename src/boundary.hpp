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

template<bool O3D> struct bdry_big {};
// LP: Can't be constexpr as std::sqrt is not constexpr without GCC extensions
template<> struct bdry_big<true>  { HOST_DEVICE static inline real value() { return RL(1.0e25); } };
template<> struct bdry_big<false> { HOST_DEVICE static inline real value() { return STD::sqrt(REAL_MAX) / RL(1.0e5); } };
#define BDRYBIG bdry_big<O3D>::value()

/*
#ifdef BHC_USE_FLOATS
#define TRIDIAG_THRESH (RL(3e-3))
#else
#define TRIDIAG_THRESH (RL(3e-6))
#endif
*/
#define TRIDIAG_THRESH (RL(1e-2))

/**
 * Copy only some of the geoacoustic data from the position-dependent array to
 * the current halfspace information. The rest of the fields are left alone from
 * the global (non-position-dependent) version.
 * [LP: FORTRAN] compiler is not accepting the copy of the whole structure at once ...
 * LP: maybe this means actually the whole struct should be copied, but he
 * only copied the elements which were needed?
 */
HOST_DEVICE inline void CopyHSInfo(HSInfo &b, const HSInfo &a)
{
    b.cP  = a.cP;
    b.cS  = a.cS;
    b.rho = a.rho;
}

/**
 * Get the top or bottom segment info (index and range interval) for range, r
 *
 * LP: t: range component of ray tangent. Endpoints of segments are handled so
 * that if the ray moves slightly along its current direction, it will remain
 * in the same segment.
 * state.lSeg: segment limits in range
 * 
 * LP: BUG: Function comment in 3D forgot to be changed when copy-pasted from 2D.
 */
template<bool O3D> HOST_DEVICE inline void GetBdrySeg(
    VEC23<O3D> x, VEC23<O3D> t, 
    BdryStateTopBot<O3D> &bds, const BdryInfoTopBot<O3D> *bdinfotb, BdryPtSmall &Bdry,
    bool isTop, bool isInit)
{
    if constexpr(O3D){
        // LP: According to the original logic, if the ray escapes the box here,
        // a warning is printed, and the other segment info is not updated. If
        // the ray escaped the box to the negative side in either dimension, or
        // the results are NaN below (should never happen), the ray is then
        // terminated in TraceRay without the current step. However, if the ray
        // escaped to the positive side, the ray is terminated as part of the
        // normal stopping conditions, including the current step.
        // On top of that weird asymmetry, of course the original logic allowed
        // for steps to the edge of the same segment they are currently in,
        // whereas we always step into the next segment according to the tangent.
        // Finally, in 2D, the altimetry / bathymetry are extended to very large
        // values so this never gets hit, whereas in Nx2D and 3D it is not so
        // almost all rays hit this case.
        // This has been changed to:
        // - An override so being at the outer edge of a segment is OK for the
        //   last segment (just for this step)
        // - If the position is actually outside the segment, print a warning
        //   and put it in the nearest valid segment.
        // - Changed the stopping condition to stop if on the boundary and
        //   pointing outwards.
        
        int32_t nx = bdinfotb->NPts.x;
        int32_t ny = bdinfotb->NPts.y;
        bds.Iseg.x = bhc::min(bhc::max(bds.Iseg.x, 0), nx-2);
        bds.Iseg.y = bhc::min(bhc::max(bds.Iseg.y, 0), ny-2);
        if(t.x >= FL(0.0)){
            while(bds.Iseg.x >= 0   && bdinfotb->bd[(bds.Iseg.x  )*ny].x.x >  x.x) --bds.Iseg.x;
            while(bds.Iseg.x >= 0   && bds.Iseg.x < nx-1 && bdinfotb->bd[(bds.Iseg.x+1)*ny].x.x <= x.x) ++bds.Iseg.x;
        }else{
            while(bds.Iseg.x < nx-1 && bdinfotb->bd[(bds.Iseg.x+1)*ny].x.x <  x.x) ++bds.Iseg.x;
            while(bds.Iseg.x >= 0   && bds.Iseg.x < nx-1 && bdinfotb->bd[(bds.Iseg.x  )*ny].x.x >= x.x) --bds.Iseg.x;
        }
        if(t.y >= FL(0.0)){
            while(bds.Iseg.y >= 0   && bdinfotb->bd[bds.Iseg.y  ].x.y >  x.y) --bds.Iseg.y;
            while(bds.Iseg.y >= 0   && bds.Iseg.y < ny-1 && bdinfotb->bd[bds.Iseg.y+1].x.y <= x.y) ++bds.Iseg.y;
        }else{
            while(bds.Iseg.y < ny-1 && bdinfotb->bd[bds.Iseg.y+1].x.y <  x.y) ++bds.Iseg.y;
            while(bds.Iseg.y >= 0   && bds.Iseg.y < ny-1 && bdinfotb->bd[bds.Iseg.y  ].x.y >= x.y) --bds.Iseg.y;
        }
        
        if(bds.Iseg.x == -1   && bdinfotb->bd[0        ].x.x == x.x) bds.Iseg.x = 0;
        if(bds.Iseg.x == nx-1 && bdinfotb->bd[(nx-1)*ny].x.x == x.x) bds.Iseg.x = nx-2;
        if(bds.Iseg.y == -1   && bdinfotb->bd[0        ].x.y == x.y) bds.Iseg.y = 0;
        if(bds.Iseg.y == ny-1 && bdinfotb->bd[ny-1     ].x.y == x.y) bds.Iseg.y = ny-2;
        
        if(bds.Iseg.x < 0 || bds.Iseg.x >= nx-1 || bds.Iseg.y < 0 || bds.Iseg.y >= ny-1){
            printf("Warning: Get%s the ray, x=(%g,%g)\n",
                isTop ? "TopSeg3D: Top altimetry undefined above" : 
                "BotSeg3D: Bottom bathymetry undefined below", x.x, x.y);
            bds.Iseg.x = bhc::min(bhc::max(bds.Iseg.x, 0), nx-2);
            bds.Iseg.y = bhc::min(bhc::max(bds.Iseg.y, 0), ny-2);
        }
        
        // segment limits in range
        bds.lSeg.x.min = bdinfotb->bd[(bds.Iseg.x  )*ny].x.x;
        bds.lSeg.x.max = bdinfotb->bd[(bds.Iseg.x+1)*ny].x.x;
        bds.lSeg.y.min = bdinfotb->bd[bds.Iseg.y  ].x.y;
        bds.lSeg.y.max = bdinfotb->bd[bds.Iseg.y+1].x.y;
        
        bds.x = bdinfotb->bd[bds.Iseg.x*ny+bds.Iseg.y].x;
        
        // printf("Iseg%s %d %d\n", isTop ? "Top" : "Bot", bds.Iseg.x+1, bds.Iseg.y+1);
        // printf("Bdryx %g,%g,%g x %g,%g,%g\n", bds.x.x, bds.x.y, bds.x.z, x.x, x.y, x.z);
        
        // identify the normal based on the active triangle of a pair
        // normal of triangle side pointing up and to the left
        vec2 tri_n = vec2(-(bds.lSeg.y.max - bds.lSeg.y.min), bds.lSeg.x.max - bds.lSeg.x.min);
        tri_n /= glm::length(tri_n);
        vec2 temp = vec2(x.x, x.y) - vec2(bds.x.x, bds.x.y);
        real over_diag_amount = glm::dot(temp, tri_n);
        // printf("temp %g,%g | tri_n %g,%g | over_diag_amount %g\n",
        //     temp.x, temp.y, tri_n.x, tri_n.y, over_diag_amount);
        if(STD::abs(over_diag_amount) > TRIDIAG_THRESH){
            bds.tridiag_pos = over_diag_amount >= RL(0.0);
        }else if(isInit){
            bds.tridiag_pos = glm::dot(vec2(t.x, t.y), tri_n) >= RL(0.0);
        }
        if(!bds.tridiag_pos){
            bds.n = bdinfotb->bd[bds.Iseg.x*ny+bds.Iseg.y].n1;
        }else{
            bds.n = bdinfotb->bd[bds.Iseg.x*ny+bds.Iseg.y].n2;
        }
        
        // if the depth is bad (a NaN) then error out
        if(!STD::isfinite(bds.x.z) || !bhc::isfinite(bds.n)){
            printf("Error: Boundary segment contains NaN!\n");
            bail();
        }
        
    }else{
        // LP: Moved from RayUpdate (TraceRay2D)
        if(!isInit &&
           !(    x.x < bds.lSeg.min || (x.x == bds.lSeg.min && t.x <  FL(0.0))
              || x.x > bds.lSeg.max || (x.x == bds.lSeg.max && t.x >= FL(0.0)) )){
            return;
        }
        
        // LP: bdinfotb->bd.x is checked for being monotonic at load time, so we can
        // linearly search out from the last position, usually only have to move
        // by 1
        int32_t n = bdinfotb->NPts;
        bds.Iseg = bhc::min(bhc::max(bds.Iseg, 0), n-2);
        if(t.x >= FL(0.0)){
            while(bds.Iseg >= 0  && bdinfotb->bd[bds.Iseg  ].x.x >  x.x) --bds.Iseg;
            while(bds.Iseg >= 0  && bds.Iseg < n-1 && bdinfotb->bd[bds.Iseg+1].x.x <= x.x) ++bds.Iseg;
        }else{
            while(bds.Iseg < n-1 && bdinfotb->bd[bds.Iseg+1].x.x <  x.x) ++bds.Iseg;
            while(bds.Iseg >= 0  && bds.Iseg < n-1 && bdinfotb->bd[bds.Iseg  ].x.x >= x.x) --bds.Iseg;
        }
        if(bds.Iseg < 0 || bds.Iseg >= n-1){
            // Iseg MUST LIE IN [0, NPts-2]
            printf("Error: Get%s the ray, r=%g\n",
                isTop ? "TopSeg: Top altimetry undefined above" : 
                "BotSeg: Bottom bathymetry undefined below", x.x);
            bail();
        }
        bds.lSeg.min = bdinfotb->bd[bds.Iseg  ].x.x;
        bds.lSeg.max = bdinfotb->bd[bds.Iseg+1].x.x;
        
        // LP: Only explicitly loaded in this function in 3D, loaded in containing
        // code in 2D
        bds.x = bdinfotb->bd[bds.Iseg].x;
        bds.n = bdinfotb->bd[bds.Iseg].n;
        
        // LP: Moved from RayInit and RayUpdate (TraceRay2D)
        if(bdinfotb->type[1] == 'L'){
            // grab the geoacoustic info for the new segment
            CopyHSInfo(Bdry.hs, bdinfotb->bd[bds.Iseg].hs);
        }
    }
}

/**
 * Does some pre-processing on the boundary points to pre-compute segment
 * lengths  (.Len),
 * tangents (.t, .nodet),
 * normals  (.n, .noden), and
 * curvatures (.kappa)
 */
template<bool O3D> inline void ComputeBdryTangentNormal(
    BdryInfoTopBot<O3D> *bd, bool isTop)
{
    typename TmplInt12<O3D>::type NPts = bd->NPts;
    vec3 tvec;
    
    if constexpr(O3D){
        
        // normals on triangle faces
        for(int32_t ix=0; ix<NPts.x - 1; ++ix){
            for(int32_t iy=0; iy<NPts.y - 1; ++iy){
                // coordinates of corner nodes, moving counter-clockwise around the rectangle
                vec3 p1 = bd->bd[(ix  )*NPts.y+iy  ].x;
                vec3 p2 = bd->bd[(ix+1)*NPts.y+iy  ].x;
                vec3 p3 = bd->bd[(ix+1)*NPts.y+iy+1].x;
                vec3 p4 = bd->bd[(ix  )*NPts.y+iy+1].x;
                
                // edges for triangle 1
                vec3 u = p2 - p1; // tangent along one edge
                vec3 v = p3 - p1; // tangent along another edge
                
                // normal vector is the cross-product of the edge tangents
                vec3 n1 = glm::cross(u, v);
                if(isTop) n1 = -n1;
                
                bd->bd[ix*NPts.y+iy].n1 = n1 / glm::length(n1); // scale to make it a unit normal
                
                // edges for triangle 2
                u = p3 - p1; // tangent along one edge
                v = p4 - p1; // tangent along another edge
                
                // normal vector is the cross-product of the edge tangents
                vec3 n2 = glm::cross(u, v);
                if(isTop) n2 = -n2;
                
                bd->bd[ix*NPts.y+iy].n2 = n2 / glm::length(n2); // scale to make it a unit normal
            }
        }
        
        // normals at nodes
        // use forward, centered, or backward difference formulas
        for(int32_t ix=0; ix<NPts.x; ++ix){
            for(int32_t iy=0; iy<NPts.y; ++iy){
                real mx, my;
                if(ix == 0){
                    mx = (bd->bd[(ix+1)*NPts.y+iy  ].x.z - bd->bd[(ix  )*NPts.y+iy  ].x.z) /
                         (bd->bd[(ix+1)*NPts.y+iy  ].x.x - bd->bd[(ix  )*NPts.y+iy  ].x.x);
                }else if(ix == NPts.x - 1){
                    mx = (bd->bd[(ix  )*NPts.y+iy  ].x.z - bd->bd[(ix-1)*NPts.y+iy  ].x.z) /
                         (bd->bd[(ix  )*NPts.y+iy  ].x.x - bd->bd[(ix-1)*NPts.y+iy  ].x.x);
                }else{
                    mx = (bd->bd[(ix+1)*NPts.y+iy  ].x.z - bd->bd[(ix-1)*NPts.y+iy  ].x.z) /
                         (bd->bd[(ix+1)*NPts.y+iy  ].x.x - bd->bd[(ix-1)*NPts.y+iy  ].x.x);
                }
                
                if(iy == 0){
                    my = (bd->bd[(ix  )*NPts.y+iy+1].x.z - bd->bd[(ix  )*NPts.y+iy  ].x.z) /
                         (bd->bd[(ix  )*NPts.y+iy+1].x.y - bd->bd[(ix  )*NPts.y+iy  ].x.y);
                }else if(iy == NPts.y - 1){
                    my = (bd->bd[(ix  )*NPts.y+iy  ].x.z - bd->bd[(ix  )*NPts.y+iy-1].x.z) /
                         (bd->bd[(ix  )*NPts.y+iy  ].x.y - bd->bd[(ix  )*NPts.y+iy-1].x.y);
                }else{
                    my = (bd->bd[(ix  )*NPts.y+iy+1].x.z - bd->bd[(ix  )*NPts.y+iy-1].x.z) /
                         (bd->bd[(ix  )*NPts.y+iy+1].x.y - bd->bd[(ix  )*NPts.y+iy-1].x.y);
                }
                
                vec3 n = vec3(-mx, -my, RL(1.0)); // this is a normal to the surface
                
                if(ix < NPts.x - 1 && iy < NPts.y - 1){
                    // xx term
                    bd->bd[(ix  )*NPts.y+iy  ].phi_xx = STD::atan2(n.z, n.x); // this is the angle at each node
                    
                    // xy term
                    tvec = bd->bd[(ix+1)*NPts.y+iy+1].x - bd->bd[(ix  )*NPts.y+iy  ].x;
                    real Len = STD::sqrt(SQ(tvec.x) + SQ(tvec.y));
                    tvec /= Len;
                    // this is the angle at each node
                    bd->bd[(ix  )*NPts.y+iy  ].phi_xy = STD::atan2(n.z, n.x * tvec.x + n.y * tvec.y);
                    
                    // yy term
                    bd->bd[(ix  )*NPts.y+iy  ].phi_yy = STD::atan2(n.z, n.y); // this is the angle at each node
                }
                
                bd->bd[(ix  )*NPts.y+iy  ].Noden_unscaled = n;
                bd->bd[(ix  )*NPts.y+iy  ].Noden = n / glm::length(n);
            }
        }
        
    }else{
        
        // LP: Moved "The boundary is also extended with a constant depth to
        // infinity to cover cases where the ray exits the domain defined by the
        // user" to ReadBoundary. The only place this is called other than there
        // is in the Init_Inline setup, which is never used and the results end
        // up the same anyway.
        
        // compute tangent and outward-pointing normal to each bottom segment
        // tBdry[0][:] = xBdry[0][1:NPts-1] - xBdry[0][0:NPts-2]
        // tBdry[1][:] = xBdry[1][1:NPts-1] - xBdry[1][0:NPts-2]
        // above caused compiler problems
        // LP: C++ obviously does not have vector slicing, but you get the idea.
        
        for(int32_t ii=0; ii<NPts-1; ++ii){
            bd->bd[ii].t  = bd->bd[ii+1].x  - bd->bd[ii].x;
            bd->bd[ii].Dx = bd->bd[ii].t[1] / bd->bd[ii].t[0]; // first derivative
            // std::cout << "Dx, t " << bd->bd[ii].Dx << " " << bd->bd[ii].x << " " 
            // << (FL(1.0) / (bd->bd[ii].x[1] / FL(500.0)) << "\n";
            
            // normalize the tangent vector
            bd->bd[ii].Len = glm::length(bd->bd[ii].t);
            bd->bd[ii].t  /= bd->bd[ii].Len;
            
            bd->bd[ii].n[0] = (isTop ? RL(1.0) : RL(-1.0)) * bd->bd[ii].t[1];
            bd->bd[ii].n[1] = (isTop ? RL(-1.0) : RL(1.0)) * bd->bd[ii].t[0];
        }
        
    }
    
    if(bd->type[0] == 'C'){
        // curvilinear option
        
        if constexpr(O3D){
            // compute derivative as centered difference between two nodes
            // compute curvatures in each segment
            
            // - sign below because the node normal = vec3(-mx, -my, RL(1.0))
            for(int32_t ix=0; ix<NPts.x; ++ix){
                for(int32_t iy=0; iy<NPts.y; ++iy){
                    real Len;
                    
                    // z_xx (difference in x of z_x)
                    bd->bd[ix*NPts.y+iy].z_xx = -(bd->bd[(ix+1)*NPts.y+iy  ].Noden_unscaled.x - bd->bd[(ix  )*NPts.y+iy  ].Noden_unscaled.x) /
                                                 (bd->bd[(ix+1)*NPts.y+iy  ].x.x              - bd->bd[(ix  )*NPts.y+iy  ].x.x);
                    
                    tvec = bd->bd[(ix+1)*NPts.y+iy  ].x - bd->bd[(ix  )*NPts.y+iy  ].x;
                    Len = STD::sqrt(SQ(tvec.x) + SQ(tvec.z));
                    // this is curvature = dphi/ds
                    bd->bd[ix*NPts.y+iy].kappa_xx = (bd->bd[(ix+1)*NPts.y+iy  ].phi_xx - bd->bd[(ix  )*NPts.y+iy  ].phi_xx) / Len;
                    
                    // z_xy (difference in y of z_x)
                    bd->bd[ix*NPts.y+iy].z_xy = -(bd->bd[(ix  )*NPts.y+iy+1].Noden_unscaled.x - bd->bd[(ix  )*NPts.y+iy  ].Noden_unscaled.x) /
                                                 (bd->bd[(ix  )*NPts.y+iy+1].x.y              - bd->bd[(ix  )*NPts.y+iy  ].x.y);
                    
                    // LP: This is overwritten by "new" below.
                    /*
                    tvec = bd->bd[(ix+1)*NPts.y+iy+1].x - bd->bd[(ix  )*NPts.y+iy  ].x;
                    Len = glm::length(tvec);
                    // this is curvature = dphi/ds
                    bd->bd[ix*NPts.y+iy].kappa_xy = (bd->bd[(ix+1)*NPts.y+iy+1].phi_xy - bd->bd[(ix  )*NPts.y+iy  ].phi_xy) / Len;
                    */
                    
                    // new
                    tvec = bd->bd[(ix  )*NPts.y+iy+1].x - bd->bd[(ix  )*NPts.y+iy  ].x;
                    Len = STD::sqrt(SQ(tvec.y) + SQ(tvec.z));
                    // this is curvature = dphi/ds
                    bd->bd[ix*NPts.y+iy].kappa_xy = (bd->bd[(ix  )*NPts.y+iy+1].phi_xx - bd->bd[(ix  )*NPts.y+iy  ].phi_xx) / Len;
                    
                    // z_yy (difference in y of z_y)
                    bd->bd[ix*NPts.y+iy].z_yy = -(bd->bd[(ix  )*NPts.y+iy+1].Noden_unscaled.y - bd->bd[(ix  )*NPts.y+iy  ].Noden_unscaled.y) /
                                                 (bd->bd[(ix  )*NPts.y+iy+1].x.y              - bd->bd[(ix  )*NPts.y+iy  ].x.y);
                    
                    tvec = bd->bd[(ix  )*NPts.y+iy+1].x - bd->bd[(ix  )*NPts.y+iy  ].x;
                    Len = STD::sqrt(SQ(tvec.y) + SQ(tvec.z));
                    // this is curvature = dphi/ds
                    bd->bd[ix*NPts.y+iy].kappa_yy = (bd->bd[(ix  )*NPts.y+iy+1].phi_xx - bd->bd[(ix  )*NPts.y+iy  ].phi_xx) / Len;
                    
                    // introduce Len factor per Eq. 4.4.18 in Cerveny's book
                    Len = glm::length(bd->bd[ix*NPts.y+iy].Noden_unscaled);
                    bd->bd[ix*NPts.y+iy].z_xx /= Len;
                    bd->bd[ix*NPts.y+iy].z_xy /= Len;
                    bd->bd[ix*NPts.y+iy].z_yy /= Len;
                }
            }
            
        }else{
        
            // compute tangent and normal at node by averaging normals on adjacent segments
            // averaging two centered differences is equivalent to forming a single centered difference of two steps ...
            for(int32_t ii=1; ii<NPts-1; ++ii){
                real sss = bd->bd[ii-1].Len / (bd->bd[ii-1].Len + bd->bd[ii].Len);
                sss = FL(0.5); // LP: BUG? Line above is overwritten.
                bd->bd[ii].Nodet = (FL(1.0) - sss) * bd->bd[ii-1].t + sss * bd->bd[ii].t;
            }
            
            bd->bd[0     ].Nodet = vec2(FL(1.0), FL(0.0)); // tangent left-end  node
            bd->bd[NPts-1].Nodet = vec2(FL(1.0), FL(0.0)); // tangent right-end node
            
            for(int32_t ii=0; ii<NPts; ++ii){
                bd->bd[ii].Noden[0] = (isTop ? RL(1.0) : RL(-1.0)) * bd->bd[ii].Nodet[1];
                bd->bd[ii].Noden[1] = (isTop ? RL(-1.0) : RL(1.0)) * bd->bd[ii].Nodet[0];
            }
        
            // compute curvature in each segment
            // LP: TODO: This allocation is not necessary, could just have two
            // variables for current and next phi. Operating on the whole array can
            // trigger compiler SIMD parallelism (AVX-512 etc.), but this is
            // unlikely to happen for atan2, and this is in one-time setup code
            // anyway.
            real *phi = allocate<real>(NPts);
            // this is the angle at each node
            for(int32_t i=0; i<NPts; ++i) phi[i] = STD::atan2(bd->bd[i].Nodet[1], bd->bd[i].Nodet[0]);
            
            for(int32_t ii=0; ii<NPts-1; ++ii){
                bd->bd[ii].kappa = (phi[ii+1] - phi[ii]) / bd->bd[ii].Len; // this is curvature = dphi/ds
                bd->bd[ii].Dxx   = (bd->bd[ii+1].Dx - bd->bd[ii].Dx) / // second derivative
                                 (bd->bd[ii+1].x[0] - bd->bd[ii].x[0]); 
                bd->bd[ii].Dss   = bd->bd[ii].Dxx * CUBE(bd->bd[ii].t[0]); // derivative in direction of tangent
                //std::cout << "kappa, Dss, Dxx " << bd->bd[ii].kappa << " " << bd->bd[ii].Dss << " " << bd->bd[ii].Dxx
                //    << " " << FL(1.0) / ((FL(8.0) / SQ(FL(1000.0))) * CUBE(STD::abs(bd->bd[ii].x[1])))
                //    << " " << bd->bd[ii].x[1] << " "
                //    << FL(-1.0) / (FL(4.0) * CUBE(bd->bd[ii].x[1]) / FL(1000000.0))
                //    << " " << bd->bd[ii].x[1] << "\n";
                
                bd->bd[ii].kappa = bd->bd[ii].Dss; // over-ride kappa !!!!!
            }
            
            deallocate(phi);
        
        }
        
    }else{
        if constexpr(O3D){
            for(int32_t ix=0; ix<NPts.x; ++ix){
                for(int32_t iy=0; iy<NPts.y; ++iy){
                    bd->bd[ix*NPts.y+iy].z_xx = RL(0.0);
                    bd->bd[ix*NPts.y+iy].z_xy = RL(0.0);
                    bd->bd[ix*NPts.y+iy].z_yy = RL(0.0);
                    
                    bd->bd[ix*NPts.y+iy].kappa_xx = RL(0.0);
                    bd->bd[ix*NPts.y+iy].kappa_xy = RL(0.0);
                    bd->bd[ix*NPts.y+iy].kappa_yy = RL(0.0);
                }
            }
        }else{
            for(int32_t i=0; i<NPts; ++i) bd->bd[i].kappa = FL(0.0);
        }
    }
}

template<bool O3D> inline void ReadBoundary(std::string FileRoot, char BdryDefMode, real BdryDepth,
    PrintFileEmu &PRTFile, BdryInfoTopBot<O3D> *bdinfotb, bool isTop,
    real freq, real fT, const AttenInfo *atten)
{
    const char *s_atibty = isTop ? "ati" : "bty";
    const char *s_ATIBTY = isTop ? "ATI" : "BTY";
    const char *s_altimetrybathymetry = isTop ? "altimetry" : "bathymetry";
    const char *s_AltimetryBathymetry = isTop ? "Altimetry" : "Bathymetry";
    const char *s_topbottom = isTop ? "top" : "bottom";
    const char *s_risesdrops = isTop ? "rises above highest" : "drops below lowest";
    
    switch(BdryDefMode){
    case '~':
    case '*':{
        if constexpr(O3D){
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
        
        LIST(BDRYFile); BDRYFile.Read(bdinfotb->type, O3D ? 1 : 2);
        if constexpr(O3D) bdinfotb->type[1] = ' ';
        switch(bdinfotb->type[0]){
        case 'R':
            if constexpr(O3D){
                PRTFile << "Regular grid for a 3D run\n";
            }else{
                PRTFile << s_atibty << "Type R not supported for 2D runs\n"; std::abort();
            }
            break;
        case 'C':
            if constexpr(O3D){
                PRTFile << "Regular grid for a 3D run (curvilinear)\n";
            }else{
                PRTFile << "Curvilinear Interpolation\n";
            }
            break;
        case 'L':
            if constexpr(O3D){
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
        
        if constexpr(O3D){
            IGNORE_UNUSED(s_risesdrops);
            
            // x values
            LIST(BDRYFile); BDRYFile.Read(bdinfotb->NPts.x);
            PRTFile << "\nNumber of " << s_altimetrybathymetry << " points in x-direction "
                << bdinfotb->NPts.x << "\n";
            
            real *Globalx = allocate<real>(std::max(bdinfotb->NPts.x, 3));
            Globalx[2] = FL(-999.9);
            LIST(BDRYFile); BDRYFile.Read(Globalx, bdinfotb->NPts.x);
            SubTab(Globalx, bdinfotb->NPts.x);
            EchoVector(Globalx, bdinfotb->NPts.x, PRTFile, Bdry_Number_to_Echo);
            // LP: BUG/TODO: This monotonic check is absent from BELLHOP3D,
            // needed for new GetBdrySeg and even implicitly for GetTop/BotSeg3D
            if(!monotonic(Globalx, bdinfotb->NPts.x)){
                std::cout << "BELLHOP:Read" << s_ATIBTY << ": " << s_AltimetryBathymetry 
                    << " X values are not monotonically increasing\n";
                std::abort();
            }
            
            // y values
            LIST(BDRYFile); BDRYFile.Read(bdinfotb->NPts.y);
            PRTFile << "\nNumber of " << s_altimetrybathymetry << " points in y-direction "
                << bdinfotb->NPts.y << "\n";
            
            real *Globaly = allocate<real>(std::max(bdinfotb->NPts.y, 3));
            Globaly[2] = FL(-999.9);
            LIST(BDRYFile); BDRYFile.Read(Globaly, bdinfotb->NPts.y);
            SubTab(Globaly, bdinfotb->NPts.y);
            EchoVector(Globaly, bdinfotb->NPts.y, PRTFile, Bdry_Number_to_Echo);
            // LP: BUG/TODO: This monotonic check is absent from BELLHOP3D,
            // needed for new GetBdrySeg and even implicitly for GetTop/BotSeg3D
            if(!monotonic(Globaly, bdinfotb->NPts.y)){
                std::cout << "BELLHOP:Read" << s_ATIBTY << ": " << s_AltimetryBathymetry 
                    << " Y values are not monotonically increasing\n";
                std::abort();
            }
            
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
                    PRTFile << "Short format (" << s_altimetrybathymetry << " only)\n";
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
                
                real sidemult = (isTop ? RL(1.0) : RL(-1.0));
                if(sidemult * bdinfotb->bd[ii].x[1] < sidemult * BdryDepth){
                    std::cout << "BELLHOP:Read" << s_ATIBTY << ": " << s_AltimetryBathymetry 
                        << " " << s_risesdrops << " point in the sound speed profile\n";
                    std::abort();
                }
            }
            
            // extend the bathymetry to +/- infinity in a piecewise constant fashion
            // LP: moved from ComputeBdryTangentNormal to be before the monotonic
            // check, as until now the first and last X were uninitialized.
            bdinfotb->bd[0               ].x[0] = -BDRYBIG;
            bdinfotb->bd[0               ].x[1] = bdinfotb->bd[1               ].x[1];
            bdinfotb->bd[0               ].hs   = bdinfotb->bd[1               ].hs;
            bdinfotb->bd[bdinfotb->NPts-1].x[0] =  BDRYBIG;
            bdinfotb->bd[bdinfotb->NPts-1].x[1] = bdinfotb->bd[bdinfotb->NPts-2].x[1];
            bdinfotb->bd[bdinfotb->NPts-1].hs   = bdinfotb->bd[bdinfotb->NPts-2].hs;
            
            // Convert ranges in km to m [LP: not the extended points]
            for(int32_t i=1; i<bdinfotb->NPts-1; ++i) bdinfotb->bd[i].x[0] *= FL(1000.0);
            
            if(!monotonic(&bdinfotb->bd[0].x.x, bdinfotb->NPts, sizeof(BdryPtFull<false>)/sizeof(real), 0)){
                std::cout << "BELLHOP:Read" << s_ATIBTY << ": " << s_AltimetryBathymetry 
                    << " ranges are not monotonically increasing\n";
                std::abort();
            }
            
        }
        
        }break;
    default:
        if constexpr(O3D){
            bdinfotb->type[0] = 'R';
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
    
    ComputeBdryTangentNormal(bdinfotb, isTop);
    
    // LP: TODO/BUG: 3D version has initialization for xTopSeg / yTopSeg here,
    // which probably also means state is carried over from one ray to the next
    
    if constexpr(!O3D){
        // convert range-dependent geoacoustic parameters from user to program units
        // LP: Moved from setup.
        if(bdinfotb->type[1] == 'L'){
            for(int32_t iSeg = 0; iSeg < bdinfotb->NPts; ++iSeg){
                // compressional wave speed
                bdinfotb->bd[iSeg].hs.cP = crci(RL(1.0e20),
                    bdinfotb->bd[iSeg].hs.alphaR, bdinfotb->bd[iSeg].hs.alphaI,
                    freq, freq, {'W', ' '}, betaPowerLaw, fT, atten, PRTFile);
                // shear         wave speed
                bdinfotb->bd[iSeg].hs.cS = crci(RL(1.0e20),
                    bdinfotb->bd[iSeg].hs.betaR, bdinfotb->bd[iSeg].hs.betaI, 
                    freq, freq, {'W', ' '}, betaPowerLaw, fT, atten, PRTFile);
            }
        }
    }else{
        IGNORE_UNUSED(freq);
        IGNORE_UNUSED(fT);
        IGNORE_UNUSED(atten);
    }
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
