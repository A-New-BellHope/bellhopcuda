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
#include "trace.hpp"
#include "jobs.hpp"

namespace bhc {

/**
 * Main ray tracing function for ray path output mode.
 */
template<bool O3D, bool R3D> HOST_DEVICE inline void MainRayMode(
    RayInitInfo &rinit, rayPt<R3D> *ray, int32_t &Nsteps, Origin<O3D, R3D> &org,
    const BdryType *ConstBdry, const BdryInfo<O3D> *bdinfo, const ReflectionInfo *refl,
    const SSPStructure *ssp, const Position *Pos, const AnglesStructure *Angles,
    const FreqInfo *freqinfo, const BeamStructure<O3D> *Beam, const BeamInfo *beaminfo)
{
    real DistBegTop, DistEndTop, DistBegBot, DistEndBot;
    SSPSegState iSeg;
    VEC23<O3D> xs, gradc;
    BdryState<O3D> bds;
    BdryType Bdry;

    if(!RayInit<O3D, R3D>(
           rinit, xs, ray[0], gradc, DistBegTop, DistBegBot, org, iSeg, bds, Bdry,
           ConstBdry, bdinfo, ssp, Pos, Angles, freqinfo, Beam, beaminfo)) {
        Nsteps = 1;
        return;
    }

    int32_t iSmallStepCtr = 0;
    int32_t is            = 0; // index for a step along the ray

    for(int32_t istep = 0; istep < MaxN - 1; ++istep) {
        is += RayUpdate<O3D, R3D>(
            ray[is], ray[is + 1], ray[is + 2], DistEndTop, DistEndBot, iSmallStepCtr, org,
            iSeg, bds, Bdry, bdinfo, refl, ssp, freqinfo, Beam, xs);
        if(RayTerminate<O3D, R3D>(
               ray[is], Nsteps, is, xs, iSmallStepCtr, DistBegTop, DistBegBot, DistEndTop,
               DistEndBot, org, bdinfo, Beam))
            break;
        if(Nsteps >= 0 && is > Nsteps) {
            Nsteps = is + 1;
            break;
        }
    }
}

template<bool O3D, bool R3D> inline void OpenRAYFile(
    LDOFile &RAYFile, std::string FileRoot, const bhcParams<O3D, R3D> &params)
{
    if(!IsRayRun(params.Beam) && !IsEigenraysRun(params.Beam)) {
        GlobalLog("OpenRAYFile not in ray trace or eigenrays mode\n");
        std::abort();
    }
    RAYFile.open(FileRoot + ".ray");
    RAYFile << params.Title << '\n';
    RAYFile << params.freqinfo->freq0 << '\n';
    RAYFile << params.Pos->NSx << params.Pos->NSy << params.Pos->NSz << '\n';
    RAYFile << params.Angles->alpha.n << params.Angles->beta.n << '\n';
    RAYFile << params.Bdry->Top.hs.Depth << '\n';
    RAYFile << params.Bdry->Bot.hs.Depth << '\n';
    RAYFile << (O3D ? "xyz" : "rz") << '\n';
}

/**
 * Compress the ray data keeping every iSkip point, points near surface or bottom, and
 * last point. Write to RAYFile.
 *
 * During an eigenray calculation, subsets of the full ray may be passed
 * These have lengths Nsteps1 vs. Nsteps for the entire ray
 *
 * The 2D version is for ray traces in (r,z) coordinates
 * The 3D version is for ray traces in (x,y,z) coordinates
 *
 * alpha0: take-off angle of this ray [LP: 2D: degrees, 3D: radians]
 */
template<bool O3D, bool R3D> void WriteRay(
    real alpha0, int32_t Nsteps1, LDOFile &RAYFile, const BdryType *Bdry,
    const Origin<O3D, R3D> &org, rayPt<R3D> *ray)
{
    // compression

    constexpr int32_t MaxNRayPoints = 500000; // this is the maximum length of the ray
                                              // vector that is written out
    int32_t n2    = 1;
    int32_t iSkip = bhc::max(Nsteps1 / MaxNRayPoints, 1);
    if constexpr(R3D) iSkip = 1; // LP: overrides line above

    for(int32_t is = 1; is < Nsteps1; ++is) {
        // ensure that we always write ray points near bdry reflections (2D only: works
        // only for flat bdry)
        if(bhc::min(
               Bdry->Bot.hs.Depth - DEP(ray[is].x), DEP(ray[is].x) - Bdry->Top.hs.Depth)
               < FL(0.2)
           || (is % iSkip) == 0 || is == Nsteps1 - 1) {
            ++n2;
            ray[n2 - 1].x = ray[is].x;
        }
    }

    // write to ray file

    if constexpr(O3D) alpha0 *= DegRad;

    RAYFile << alpha0 << '\n';
    RAYFile << n2 << ray[Nsteps1 - 1].NumTopBnc << ray[Nsteps1 - 1].NumBotBnc << '\n';

    for(int32_t is = 0; is < n2; ++is) { RAYFile << RayToOceanX(ray[is].x, org) << '\n'; }
}

template<bool O3D, bool R3D> inline void InitRayMode(
    RayInfo<O3D, R3D> *rayinfo, const bhcParams<O3D, R3D> &params)
{
    rayinfo->NRays     = GetNumJobs<O3D>(params.Pos, params.Angles);
    rayinfo->MaxPoints = bhc::min((uint32_t)MaxN * (uint32_t)rayinfo->NRays, 100000000u);
    rayinfo->NPoints   = 0;
    checkallocate(rayinfo->raymem, rayinfo->MaxPoints);
    checkallocate(rayinfo->results, rayinfo->NRays);
    memset(rayinfo->results, 0, rayinfo->NRays * sizeof(RayResult<O3D, R3D>)); // Clear
                                                                               // because
                                                                               // will
                                                                               // check
                                                                               // pointers
}

template<bool O3D, bool R3D> inline void FinalizeRayMode(
    RayInfo<O3D, R3D> *rayinfo, std::string FileRoot, const bhcParams<O3D, R3D> &params)
{
    LDOFile RAYFile;
    OpenRAYFile<O3D, R3D>(RAYFile, FileRoot, params);
    for(int r = 0; r < rayinfo->NRays; ++r) {
        RayResult<O3D, R3D> *res = &rayinfo->results[r];
        if(res->ray == nullptr) continue;
        WriteRay<O3D, R3D>(
            res->SrcDeclAngle, res->Nsteps, RAYFile, params.Bdry, res->org, res->ray);
    }
}

template<bool O3D, bool R3D> inline bool IsRayCopyMode(const RayInfo<O3D, R3D> *rayinfo)
{
    return (size_t)rayinfo->MaxPoints < (size_t)MaxN * (size_t)rayinfo->NRays;
}

template<bool O3D, bool R3D> inline bool RunRay(
    RayInfo<O3D, R3D> *rayinfo, const bhcParams<O3D, R3D> &params, rayPt<R3D> *localmem,
    int32_t job, RayInitInfo &rinit, int32_t &Nsteps)
{
    rayPt<R3D> *ray;
    if(IsRayCopyMode<O3D, R3D>(rayinfo)) {
        ray = localmem;
    } else {
        ray = &rayinfo->raymem[job * MaxN];
    }
    memset(ray, 0xFE, MaxN * sizeof(rayPt<R3D>)); // Set to garbage values for debugging

    Origin<O3D, R3D> org;
    MainRayMode<O3D, R3D>(
        rinit, ray, Nsteps, org, params.Bdry, params.bdinfo, params.refl, params.ssp,
        params.Pos, params.Angles, params.freqinfo, params.Beam, params.beaminfo);

    bool ret = true;
    if(IsRayCopyMode<O3D, R3D>(rayinfo)) {
        uint32_t p = AtomicFetchAdd(&rayinfo->NPoints, (uint32_t)Nsteps);
        if(p + Nsteps > rayinfo->MaxPoints) {
            GlobalLog("Ran out of memory for rays\n");
            rayinfo->results[job].ray = nullptr;
            ret                       = false;
        } else {
            rayinfo->results[job].ray = &rayinfo->raymem[p];
            memcpy(rayinfo->results[job].ray, localmem, Nsteps * sizeof(rayPt<R3D>));
        }
    } else {
        rayinfo->results[job].ray = ray;
    }
    rayinfo->results[job].org          = org;
    rayinfo->results[job].SrcDeclAngle = rinit.SrcDeclAngle;
    rayinfo->results[job].Nsteps       = Nsteps;

    return ret;
}

} // namespace bhc
