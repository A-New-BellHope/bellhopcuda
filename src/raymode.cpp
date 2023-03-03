/*
bellhopcxx / bellhopcuda - C++/CUDA port of BELLHOP underwater acoustics simulator
Copyright (C) 2021-2023 The Regents of the University of California
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
#include "raymode.hpp"

namespace bhc {

template<bool O3D, bool R3D> inline void OpenRAYFile(
    LDOFile &RAYFile, std::string FileRoot, const bhcParams<O3D, R3D> &params)
{
    if(!IsRayRun(params.Beam) && !IsEigenraysRun(params.Beam)) {
        EXTERR("OpenRAYFile not in ray trace or eigenrays mode");
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
 * last point.
 *
 * The 2D version is for ray traces in (r,z) coordinates
 * The 3D version is for ray traces in (x,y,z) coordinates
 */
template<bool O3D, bool R3D> inline void CompressRay(
    RayResult<O3D, R3D> *res, const BdryType *Bdry)
{
    // this is the maximum length of the ray vector that is written out
    constexpr int32_t MaxNRayPoints = 500000;

    // compression
    // LP: This is silly for two reasons:
    // 1) MaxN (maximum number of steps for a ray) is 100000, but MaxNRayPoints
    //    is 500000. Therefore iSkip will always be 1, and the whole vector will
    //    always be written.
    // 2) Even if these constants were changed, the formula for iSkip is not
    //    ideal: iSkip will only become 2 once the number of steps in the ray is
    //    more than 2x MaxNRayPoints. If it's less than this, it'll just be
    //    truncated, which is arguably worse than skipping every other step.
    // So we'll just make sure this doesn't run unless the constants are changed.
    if constexpr(MaxN > MaxNRayPoints) {
        int32_t n2    = 1;
        int32_t iSkip = bhc::max(res->Nsteps / MaxNRayPoints, 1);
        if constexpr(R3D) iSkip = 1; // LP: overrides line above

        for(int32_t is = 1; is < res->Nsteps; ++is) {
            // ensure that we always write ray points near bdry reflections (2D only:
            // works only for flat bdry)
            if(bhc::min(
                   Bdry->Bot.hs.Depth - DEP(res->ray[is].x),
                   DEP(res->ray[is].x) - Bdry->Top.hs.Depth)
                   < FL(0.2)
               || (is % iSkip) == 0 || is == res->Nsteps - 1) {
                ++n2;
                res->ray[n2 - 1].x = res->ray[is].x;
            }
        }
        res->Nsteps = n2;
    }
}

/**
 * Write to RAYFile.
 */
template<bool O3D, bool R3D> inline void WriteRay(
    LDOFile &RAYFile, const RayResult<O3D, R3D> *res)
{
    // take-off angle of this ray [LP: 2D: degrees, 3D: radians]
    real alpha0 = res->SrcDeclAngle;
    if constexpr(O3D) alpha0 *= DegRad;
    RAYFile << alpha0 << '\n';

    RAYFile << res->Nsteps;
    RAYFile << res->ray[res->Nsteps - 1].NumTopBnc;
    RAYFile << res->ray[res->Nsteps - 1].NumBotBnc << '\n';
    for(int32_t is = 0; is < res->Nsteps; ++is) {
        RAYFile << RayToOceanX(res->ray[is].x, res->org) << '\n';
    }
}

/**
 * neigen is 0 when not in eigenrays mode (i.e. in ray trace mode).
 */
template<bool O3D, bool R3D> void InitRayMode(
    RayInfo<O3D, R3D> *rayinfo, const bhcParams<O3D, R3D> &params, uint32_t neigen)
{
#warning TODO InitRayMode
    rayinfo->NRays     = neigen > 0 ? neigen : GetNumJobs<O3D>(params.Pos, params.Angles);
    rayinfo->MaxPoints = bhc::min((uint32_t)MaxN * (uint32_t)rayinfo->NRays, 100000000u);
    rayinfo->NPoints   = 0;
    trackallocate(params, "rays", rayinfo->raymem, rayinfo->MaxPoints);
    trackallocate(params, "rays", rayinfo->results, rayinfo->NRays);
    // Clear because will check pointers
    memset(rayinfo->results, 0, rayinfo->NRays * sizeof(RayResult<O3D, R3D>));
}

#if BHC_ENABLE_2D
template void InitRayMode<false, false>(
    RayInfo<false, false> *rayinfo, const bhcParams<false, false> &params,
    uint32_t neigen);
#endif
#if BHC_ENABLE_NX2D
template void InitRayMode<true, false>(
    RayInfo<true, false> *rayinfo, const bhcParams<true, false> &params, uint32_t neigen);
#endif
#if BHC_ENABLE_3D
template void InitRayMode<true, true>(
    RayInfo<true, true> *rayinfo, const bhcParams<true, true> &params, uint32_t neigen);
#endif

template<bool O3D, bool R3D> void PostProcessRays(
    const bhcParams<O3D, R3D> &params, RayInfo<O3D, R3D> *rayinfo)
{
    for(int r = 0; r < rayinfo->NRays; ++r) {
        RayResult<O3D, R3D> *res = &rayinfo->results[r];
        if(res->ray == nullptr) continue;
        CompressRay<O3D, R3D>(res, params.Bdry);
    }
}

#if BHC_ENABLE_2D
template void PostProcessRays<false, false>(
    const bhcParams<false, false> &params, RayInfo<false, false> *rayinfo);
#endif
#if BHC_ENABLE_NX2D
template void PostProcessRays<true, false>(
    const bhcParams<true, false> &params, RayInfo<true, false> *rayinfo);
#endif
#if BHC_ENABLE_3D
template void PostProcessRays<true, true>(
    const bhcParams<true, true> &params, RayInfo<true, true> *rayinfo);
#endif

template<bool O3D, bool R3D> void WriteOutRays(
    const bhcParams<O3D, R3D> &params, const RayInfo<O3D, R3D> *rayinfo)
{
    LDOFile RAYFile;
    OpenRAYFile<O3D, R3D>(RAYFile, GetInternal(params)->FileRoot, params);
    for(int r = 0; r < rayinfo->NRays; ++r) {
        const RayResult<O3D, R3D> *res = &rayinfo->results[r];
        if(res->ray == nullptr) continue;
        WriteRay<O3D, R3D>(RAYFile, res);
    }
}

#if BHC_ENABLE_2D
template void WriteOutRays<false, false>(
    const bhcParams<false, false> &params, const RayInfo<false, false> *rayinfo);
#endif
#if BHC_ENABLE_NX2D
template void WriteOutRays<true, false>(
    const bhcParams<true, false> &params, const RayInfo<true, false> *rayinfo);
#endif
#if BHC_ENABLE_3D
template void WriteOutRays<true, true>(
    const bhcParams<true, true> &params, const RayInfo<true, true> *rayinfo);
#endif

} // namespace bhc
