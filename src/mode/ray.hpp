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
#include "../common_setup.hpp"
#include "modemodule.hpp"

namespace bhc { namespace mode {

template<bool O3D, bool R3D> bool RunRay(
    RayInfo<O3D, R3D> *rayinfo, const bhcParams<O3D> &params, int32_t job, int32_t worker,
    RayInitInfo &rinit, int32_t &Nsteps, ErrState *errState);
extern template bool RunRay<false, false>(
    RayInfo<false, false> *rayinfo, const bhcParams<false> &params, int32_t job,
    int32_t worker, RayInitInfo &rinit, int32_t &Nsteps, ErrState *errState);
extern template bool RunRay<true, false>(
    RayInfo<true, false> *rayinfo, const bhcParams<true> &params, int32_t job,
    int32_t worker, RayInitInfo &rinit, int32_t &Nsteps, ErrState *errState);
extern template bool RunRay<true, true>(
    RayInfo<true, true> *rayinfo, const bhcParams<true> &params, int32_t job,
    int32_t worker, RayInitInfo &rinit, int32_t &Nsteps, ErrState *errState);

template<bool O3D, bool R3D> void RunRayMode(
    bhcParams<O3D> &params, bhcOutputs<O3D, R3D> &outputs);
extern template void RunRayMode<false, false>(
    bhcParams<false> &params, bhcOutputs<false, false> &outputs);
extern template void RunRayMode<true, false>(
    bhcParams<true> &params, bhcOutputs<true, false> &outputs);
extern template void RunRayMode<true, true>(
    bhcParams<true> &params, bhcOutputs<true, true> &outputs);

template<bool O3D, bool R3D> void ReadOutRay(
    bhcParams<O3D> &params, bhcOutputs<O3D, R3D> &outputs, const char *FileRoot);
extern template void ReadOutRay<false, false>(
    bhcParams<false> &params, bhcOutputs<false, false> &outputs, const char *FileRoot);
extern template void ReadOutRay<true, false>(
    bhcParams<true> &params, bhcOutputs<true, false> &outputs, const char *FileRoot);
extern template void ReadOutRay<true, true>(
    bhcParams<true> &params, bhcOutputs<true, true> &outputs, const char *FileRoot);

template<bool O3D, bool R3D> class Ray : public ModeModule<O3D, R3D> {
public:
    Ray() {}
    virtual ~Ray() {}

    virtual void Init(bhcOutputs<O3D, R3D> &outputs) const override
    {
        outputs.rayinfo->results    = nullptr;
        outputs.rayinfo->RayMem     = nullptr;
        outputs.rayinfo->WorkRayMem = nullptr;

        outputs.rayinfo->RayMemCapacity  = 0;
        outputs.rayinfo->RayMemPoints    = 0;
        outputs.rayinfo->MaxPointsPerRay = 0;
        outputs.rayinfo->NRays           = 0;
    }

    virtual void Preprocess(
        bhcParams<O3D> &params, bhcOutputs<O3D, R3D> &outputs) const override
    {
        RayInfo<O3D, R3D> *rayinfo = outputs.rayinfo;

        trackdeallocate(params, rayinfo->RayMem);
        trackdeallocate(params, rayinfo->WorkRayMem);
        rayinfo->NRays = IsEigenraysRun(params.Beam)
            ? outputs.eigen->neigen
            : GetNumJobs<O3D>(params.Pos, params.Angles);
        trackallocate(params, "ray metadata", rayinfo->results, rayinfo->NRays);
        // Clear because will check pointers
        memset(rayinfo->results, 0, rayinfo->NRays * sizeof(RayResult<O3D, R3D>));

        rayinfo->MaxPointsPerRay = MaxN;
        rayinfo->isCopyMode      = false;
        size_t needtotalsize = (size_t)rayinfo->NRays * (size_t)MaxN * sizeof(rayPt<R3D>);
        if(GetInternal(params)->usedMemory + needtotalsize
           <= GetInternal(params)->maxMemory) {
            rayinfo->RayMemCapacity = (size_t)rayinfo->NRays * (size_t)MaxN;
        } else if(GetInternal(params)->useRayCopyMode) {
            trackallocate(
                params, "work rays for copy mode", rayinfo->WorkRayMem,
                GetInternal(params)->numThreads * MaxN);
            rayinfo->RayMemCapacity = (GetInternal(params)->maxMemory
                                       - GetInternal(params)->usedMemory)
                / sizeof(rayPt<R3D>);
            rayinfo->isCopyMode = true;
        } else {
            rayinfo->MaxPointsPerRay = (int32_t)std::min(
                (GetInternal(params)->maxMemory - GetInternal(params)->usedMemory)
                    / ((size_t)rayinfo->NRays * sizeof(rayPt<R3D>)),
                (size_t)0x7FFFFFFF);
            if(rayinfo->MaxPointsPerRay == 0) {
                EXTERR("Insufficient memory to allocate any rays at all");
            } else if(rayinfo->MaxPointsPerRay < 500) {
                EXTWARN(
                    "There is only enough memory to allocate %d points per ray",
                    rayinfo->MaxPointsPerRay);
            }
            rayinfo->RayMemCapacity = (size_t)rayinfo->NRays
                * (size_t)rayinfo->MaxPointsPerRay;
        }
        trackallocate(params, "rays", rayinfo->RayMem, rayinfo->RayMemCapacity);
        rayinfo->RayMemPoints = 0;
    }

    virtual void Run(bhcParams<O3D> &params, bhcOutputs<O3D, R3D> &outputs) const override
    {
        RunRayMode<O3D, R3D>(params, outputs);
    }

    virtual void Postprocess(
        bhcParams<O3D> &params, bhcOutputs<O3D, R3D> &outputs) const override
    {
        RayInfo<O3D, R3D> *rayinfo = outputs.rayinfo;
        for(int r = 0; r < rayinfo->NRays; ++r) {
            RayResult<O3D, R3D> *res = &rayinfo->results[r];
            if(res->ray == nullptr) continue;
            CompressRay(res, params.Bdry);
        }
    }

    virtual void Writeout(
        const bhcParams<O3D> &params, const bhcOutputs<O3D, R3D> &outputs) const override
    {
        RayInfo<O3D, R3D> *rayinfo = outputs.rayinfo;
        LDOFile RAYFile;
        OpenRAYFile(RAYFile, GetInternal(params)->FileRoot, params);
        for(int r = 0; r < rayinfo->NRays; ++r) {
            const RayResult<O3D, R3D> *res = &rayinfo->results[r];
            if(res->ray == nullptr) continue;
            WriteRay(RAYFile, res);
        }
    }

    virtual void Readout(
        bhcParams<O3D> &params, bhcOutputs<O3D, R3D> &outputs,
        const char *FileRoot) const override
    {
        ReadOutRay<O3D, R3D>(params, outputs, FileRoot);
    }

    virtual void Finalize(
        bhcParams<O3D> &params, bhcOutputs<O3D, R3D> &outputs) const override
    {
        trackdeallocate(params, outputs.rayinfo->results);
        trackdeallocate(params, outputs.rayinfo->RayMem);
        trackdeallocate(params, outputs.rayinfo->WorkRayMem);
    }

private:
    // LP: These are small enough that it's not really necessary to compile
    // them separately.

    inline void OpenRAYFile(
        LDOFile &RAYFile, std::string FileRoot, const bhcParams<O3D> &params) const
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
    inline void CompressRay(RayResult<O3D, R3D> *res, const BdryType *Bdry) const
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
    inline void WriteRay(LDOFile &RAYFile, const RayResult<O3D, R3D> *res) const
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
};

}} // namespace bhc::mode
