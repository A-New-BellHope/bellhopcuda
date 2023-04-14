/*
bellhopcxx / bellhopcuda - C++/CUDA port of BELLHOP(3D) underwater acoustics simulator
Copyright (C) 2021-2023 The Regents of the University of California
Marine Physical Lab at Scripps Oceanography, c/o Jules Jaffe, jjaffe@ucsd.edu
Based on BELLHOP / BELLHOP3D, which is Copyright (C) 1983-2022 Michael B. Porter

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
#include "common_run.hpp"

namespace bhc {

constexpr int32_t Bdry_Number_to_Echo = 21;

template<bool O3D> struct bdry_big {};
// LP: Can't be constexpr as std::sqrt is not constexpr without GCC extensions
template<> struct bdry_big<true> {
    HOST_DEVICE static inline real value() { return RL(1.0e25); }
};
template<> struct bdry_big<false> {
    HOST_DEVICE static inline real value() { return STD::sqrt(REAL_MAX) / RL(1.0e5); }
};
#define BDRYBIG bdry_big<O3D>::value()

#ifdef BHC_USE_FLOATS
#define TRIDIAG_THRESH (RL(3e-3))
#else
#define TRIDIAG_THRESH (RL(3e-6))
#endif

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
 * Get the top or bottom segment info (index and range interval) for range, r,
 * or XY position, x
 *
 * LP: t: range component of ray tangent. Endpoints of segments are handled so
 * that if the ray moves slightly along its current direction, it will remain
 * in the same segment.
 * state.lSeg: segment limits in range
 */
template<bool O3D> HOST_DEVICE inline void GetBdrySeg(
    VEC23<O3D> x, VEC23<O3D> t, BdryStateTopBot<O3D> &bds,
    const BdryInfoTopBot<O3D> *bdinfotb, BdryPtSmall &Bdry, bool isTop, bool isInit,
    ErrState *errState)
{
    if constexpr(O3D) {
        // LP: See discussion of changes in Fortran version readme.

        int32_t nx = bdinfotb->NPts.x;
        int32_t ny = bdinfotb->NPts.y;
        bds.Iseg.x = bhc::min(bhc::max(bds.Iseg.x, 0), nx - 2);
        bds.Iseg.y = bhc::min(bhc::max(bds.Iseg.y, 0), ny - 2);
        if(t.x >= FL(0.0)) {
            while(bds.Iseg.x >= 0 && bdinfotb->bd[(bds.Iseg.x) * ny].x.x > x.x)
                --bds.Iseg.x;
            while(bds.Iseg.x >= 0 && bds.Iseg.x < nx - 1
                  && bdinfotb->bd[(bds.Iseg.x + 1) * ny].x.x <= x.x)
                ++bds.Iseg.x;
        } else {
            while(bds.Iseg.x < nx - 1 && bdinfotb->bd[(bds.Iseg.x + 1) * ny].x.x < x.x)
                ++bds.Iseg.x;
            while(bds.Iseg.x >= 0 && bds.Iseg.x < nx - 1
                  && bdinfotb->bd[(bds.Iseg.x) * ny].x.x >= x.x)
                --bds.Iseg.x;
        }
        if(t.y >= FL(0.0)) {
            while(bds.Iseg.y >= 0 && bdinfotb->bd[bds.Iseg.y].x.y > x.y) --bds.Iseg.y;
            while(bds.Iseg.y >= 0 && bds.Iseg.y < ny - 1
                  && bdinfotb->bd[bds.Iseg.y + 1].x.y <= x.y)
                ++bds.Iseg.y;
        } else {
            while(bds.Iseg.y < ny - 1 && bdinfotb->bd[bds.Iseg.y + 1].x.y < x.y)
                ++bds.Iseg.y;
            while(bds.Iseg.y >= 0 && bds.Iseg.y < ny - 1
                  && bdinfotb->bd[bds.Iseg.y].x.y >= x.y)
                --bds.Iseg.y;
        }

        if(bds.Iseg.x == -1 && bdinfotb->bd[0].x.x == x.x) bds.Iseg.x = 0;
        if(bds.Iseg.x == nx - 1 && bdinfotb->bd[(nx - 1) * ny].x.x == x.x)
            bds.Iseg.x = nx - 2;
        if(bds.Iseg.y == -1 && bdinfotb->bd[0].x.y == x.y) bds.Iseg.y = 0;
        if(bds.Iseg.y == ny - 1 && bdinfotb->bd[ny - 1].x.y == x.y) bds.Iseg.y = ny - 2;

        if(bds.Iseg.x < 0 || bds.Iseg.x >= nx - 1 || bds.Iseg.y < 0
           || bds.Iseg.y >= ny - 1) {
            RunWarning(
                errState,
                isTop ? BHC_WARN_OUTSIDE_ALTIMETRY : BHC_WARN_OUTSIDE_BATHYMETRY);
            /*
            printf(
                "Warning: Get%s the ray, x=(%g,%g)\n",
                isTop ? "TopSeg3D: Top altimetry undefined above"
                      : "BotSeg3D: Bottom bathymetry undefined below",
                x.x, x.y);
            */
            bds.Iseg.x = bhc::min(bhc::max(bds.Iseg.x, 0), nx - 2);
            bds.Iseg.y = bhc::min(bhc::max(bds.Iseg.y, 0), ny - 2);
        }

        // segment limits in range
        bds.lSeg.x.min = bdinfotb->bd[(bds.Iseg.x) * ny].x.x;
        bds.lSeg.x.max = bdinfotb->bd[(bds.Iseg.x + 1) * ny].x.x;
        bds.lSeg.y.min = bdinfotb->bd[bds.Iseg.y].x.y;
        bds.lSeg.y.max = bdinfotb->bd[bds.Iseg.y + 1].x.y;

        bds.x    = bdinfotb->bd[bds.Iseg.x * ny + bds.Iseg.y].x;
        bds.xmid = (bds.x + bdinfotb->bd[(bds.Iseg.x + 1) * ny + (bds.Iseg.y + 1)].x)
            * RL(0.5);

        // printf("Iseg%s %d %d\n", isTop ? "Top" : "Bot", bds.Iseg.x+1, bds.Iseg.y+1);
        // printf("Bdryx %g,%g,%g x %g,%g,%g\n", bds.x.x, bds.x.y, bds.x.z, x.x, x.y,
        // x.z);

        // identify the normal based on the active triangle of a pair
        // normal of triangle side pointing up and to the left
        vec2 tri_n
            = vec2(-(bds.lSeg.y.max - bds.lSeg.y.min), bds.lSeg.x.max - bds.lSeg.x.min);
        tri_n /= glm::length(tri_n);
        vec2 temp             = vec2(x.x, x.y) - vec2(bds.xmid.x, bds.xmid.y);
        real over_diag_amount = glm::dot(temp, tri_n);
        bds.td.onEdge         = STD::abs(over_diag_amount) < TRIDIAG_THRESH;
        // printf("temp %g,%g | tri_n %g,%g | over_diag_amount %g\n",
        //     temp.x, temp.y, tri_n.x, tri_n.y, over_diag_amount);
        if(!isInit && bds.td.justSteppedTo) {
            bds.td.side = bds.td.outgoingSide;
        } else if(!isInit && bds.td.onEdge) {
            bds.td.side = glm::dot(XYCOMP(t), tri_n) >= RL(0.0);
        } else {
            bds.td.side = over_diag_amount >= RL(0.0);
        }
        bds.td.justSteppedTo = false;
        if(!bds.td.side) {
            bds.n = bdinfotb->bd[bds.Iseg.x * ny + bds.Iseg.y].n1;
        } else {
            bds.n = bdinfotb->bd[bds.Iseg.x * ny + bds.Iseg.y].n2;
        }

        // if the depth is bad (a NaN) then error out
        if(!STD::isfinite(bds.x.z) || !bhc::isfinite(bds.n)) {
            RunError(errState, BHC_ERR_BOUNDARY_SEG_CONTAINS_NAN);
        }

    } else {
        // LP: Moved from RayUpdate (TraceRay2D)
        if(!isInit
           && !(
               x.x < bds.lSeg.min || (x.x == bds.lSeg.min && t.x < FL(0.0))
               || x.x > bds.lSeg.max || (x.x == bds.lSeg.max && t.x >= FL(0.0)))) {
            return;
        }

        // LP: bdinfotb->bd.x is checked for being monotonic at load time, so we can
        // linearly search out from the last position, usually only have to move
        // by 1
        int32_t n = bdinfotb->NPts;
        bds.Iseg  = bhc::min(bhc::max(bds.Iseg, 0), n - 2);
        if(t.x >= FL(0.0)) {
            while(bds.Iseg >= 0 && bdinfotb->bd[bds.Iseg].x.x > x.x) --bds.Iseg;
            while(bds.Iseg >= 0 && bds.Iseg < n - 1
                  && bdinfotb->bd[bds.Iseg + 1].x.x <= x.x)
                ++bds.Iseg;
        } else {
            while(bds.Iseg < n - 1 && bdinfotb->bd[bds.Iseg + 1].x.x < x.x) ++bds.Iseg;
            while(bds.Iseg >= 0 && bds.Iseg < n - 1 && bdinfotb->bd[bds.Iseg].x.x >= x.x)
                --bds.Iseg;
        }
        if(bds.Iseg < 0 || bds.Iseg >= n - 1) {
            // Iseg MUST LIE IN [0, NPts-2]
            RunError(
                errState, isTop ? BHC_ERR_OUTSIDE_ALTIMETRY : BHC_ERR_OUTSIDE_BATHYMETRY);
            /*
            printf(
                "Error: Get%s the ray, r=%g\n",
                isTop ? "TopSeg: Top altimetry undefined above"
                      : "BotSeg: Bottom bathymetry undefined below",
                x.x);
            */
            bds.Iseg = 0;
        }
        bds.lSeg.min = bdinfotb->bd[bds.Iseg].x.x;
        bds.lSeg.max = bdinfotb->bd[bds.Iseg + 1].x.x;

        // LP: Only explicitly loaded in this function in 3D, loaded in containing
        // code in 2D
        bds.x = bdinfotb->bd[bds.Iseg].x;
        // bds.xmid = (bds.x + bdinfotb->bd[bds.Iseg+1].x) * RL(0.5);
        bds.n = bdinfotb->bd[bds.Iseg].n;

        // LP: Moved from RayInit and RayUpdate (TraceRay2D)
        if(bdinfotb->type[1] == 'L') {
            // grab the geoacoustic info for the new segment
            CopyHSInfo(Bdry.hs, bdinfotb->bd[bds.Iseg].hs);
        }
    }
}

} // namespace bhc
