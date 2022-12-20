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
#include "raymode.hpp"
#include "tlmode.hpp"
#include "eigenrays.hpp"
#include "arrivals.hpp"
#include "ssp.hpp"

#include <atomic>
#include <mutex>

namespace bhc {

extern std::atomic<int32_t> sharedJobID;
extern std::mutex exceptionMutex;
extern std::string exceptionStr;

#ifdef BHC_BUILD_CUDA
extern int m_gpu, d_warp, d_maxthreads, d_multiprocs;
#endif

template<bool O3D, bool R3D> inline void InitSelectedMode(
    bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs, bool singlethread)
{
    // this is always called from run_* so update intermediate params
    PrintFileEmu &PRTFile = *(PrintFileEmu *)params.internal;
    UpdateSSP(params.freqinfo->freq0, params.fT, params.ssp, PRTFile, params.atten, O3D);

    // Common
    int32_t ns = params.Pos->NSx * params.Pos->NSy * params.Pos->NSz;
    BASSERT(ns >= 1);

    // irregular or rectilinear grid
    params.Pos->NRz_per_range = IsIrregularGrid(params.Beam) ? 1 : params.Pos->NRz;

    // Mode specific
    if(IsRayRun(params.Beam)) {
        InitRayMode<O3D, R3D>(outputs.rayinfo, params, 0);
    } else if(IsTLRun(params.Beam)) {
        InitTLMode(outputs.uAllSources, params.Pos);
    } else if(IsEigenraysRun(params.Beam)) {
        InitEigenMode(outputs.eigen);
    } else if(IsArrivalsRun(params.Beam)) {
        InitArrivalsMode(
            outputs.arrinfo, singlethread, O3D, params.Pos,
            *(PrintFileEmu *)params.internal);
    } else {
        GlobalLog("Invalid RunType %c\n", params.Beam->RunType[0]);
        std::abort();
    }

    if(!IsRayRun(params.Beam)) { PreRun_Influence<O3D, R3D>(params.Beam, params.Pos); }
}

template<typename CFG, bool O3D, bool R3D> void FieldModesWorker(
    bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs);

template<typename CFG, bool O3D, bool R3D> void RunFieldModesImpl(
    bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs, uint32_t cores);

template<char RT, char IT, bool O3D, bool R3D> inline void RunFieldModesSelSSP(
    bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs, uint32_t cores)
{
    char st = params.ssp->Type;
    if(st == 'N') {
#ifdef BHC_SSP_ENABLE_N2LINEAR
        RunFieldModesImpl<CfgSel<RT, IT, 'N'>, O3D, R3D>(params, outputs, cores);
#else
        GlobalLog("N2-linear SSP (ssp->Type == 'N') was not enabled at compile time!");
        bail();
#endif
    } else if(st == 'C') {
#ifdef BHC_SSP_ENABLE_CLINEAR
        RunFieldModesImpl<CfgSel<RT, IT, 'C'>, O3D, R3D>(params, outputs, cores);
#else
        GlobalLog("C-linear SSP (ssp->Type == 'C') was not enabled at compile time!");
        bail();
#endif
    } else if(st == 'S') {
#ifdef BHC_SSP_ENABLE_CUBIC
        RunFieldModesImpl<CfgSel<RT, IT, 'S'>, O3D, R3D>(params, outputs, cores);
#else
        GlobalLog("Cubic spline SSP (ssp->Type == 'S') was not enabled at compile time!");
        bail();
#endif
    } else if(st == 'P') {
#ifdef BHC_SSP_ENABLE_PCHIP
#ifdef BHC_LIMIT_FEATURES
        if constexpr(!O3D) {
#endif
            RunFieldModesImpl<CfgSel<RT, IT, 'P'>, O3D, R3D>(params, outputs, cores);
#ifdef BHC_LIMIT_FEATURES
        } else {
            GlobalLog("Nx2D or 3D PCHIP SSP not supported"
                      "because BHC_LIMIT_FEATURES enabled!");
            bail();
        }
#endif
#else
        GlobalLog("PCHIP SSP (ssp->Type == 'P') was not enabled at compile time!");
        bail();
#endif
    } else if(st == 'Q') {
#ifdef BHC_SSP_ENABLE_QUAD
        if constexpr(!O3D) {
            RunFieldModesImpl<CfgSel<RT, IT, 'Q'>, O3D, R3D>(params, outputs, cores);
        } else {
            GlobalLog("Quad SSP not supported in Nx2D or 3D mode!");
            bail();
        }
#else
        GlobalLog("Quad SSP (ssp->Type == 'Q') was not enabled at compile time!");
        bail();
#endif
    } else if(st == 'H') {
#ifdef BHC_SSP_ENABLE_HEXAHEDRAL
        if constexpr(O3D) {
            RunFieldModesImpl<CfgSel<RT, IT, 'H'>, O3D, R3D>(params, outputs, cores);
        } else {
            GlobalLog("Hexahedral SSP not supported in 2D mode!");
            bail();
        }
#else
        GlobalLog("Hexahedral SSP (ssp->Type == 'H') was not enabled at compile time!");
        bail();
#endif
    } else if(st == 'A') {
#ifdef BHC_SSP_ENABLE_ANALYTIC
        RunFieldModesImpl<CfgSel<RT, IT, 'A'>, O3D, R3D>(params, outputs, cores);
#else
        GlobalLog("Analytic SSP (ssp->Type == 'A') was not enabled at compile time!");
        bail();
#endif
    } else {
        GlobalLog("Invalid ssp->Type %c!", st);
        bail();
    }
}

template<char IT, bool O3D, bool R3D> inline void RunFieldModesSelRun(
    bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs, uint32_t cores)
{
    char rt = params.Beam->RunType[0];
    if(rt == 'C' || rt == 'S' || rt == 'I') {
#ifdef BHC_RUN_ENABLE_TL
        RunFieldModesSelSSP<'C', IT, O3D, R3D>(params, outputs, cores);
#else
        GlobalLog("Transmission loss runs (Beam->RunType[0] == 'C', 'S', or 'I') "
                  "were not enabled at compile time!");
        bail();
#endif
    } else if(rt == 'E') {
#ifdef BHC_RUN_ENABLE_EIGENRAYS
        if constexpr(InflType<IT>::IsCerveny()) {
            GlobalLog("Cerveny influence does not support eigenrays!");
            bail();
        } else {
            RunFieldModesSelSSP<'E', IT, O3D, R3D>(params, outputs, cores);
        }
#else
        GlobalLog("Eigenrays runs (Beam->RunType[0] == 'E') "
                  "were not enabled at compile time!");
        bail();
#endif
    } else if(rt == 'A' || rt == 'a') {
#ifdef BHC_RUN_ENABLE_ARRIVALS
        if constexpr(InflType<IT>::IsCerveny()) {
            GlobalLog("Cerveny influence does not support arrivals!");
            bail();
        } else {
            RunFieldModesSelSSP<'A', IT, O3D, R3D>(params, outputs, cores);
        }
#else
        GlobalLog("Arrivals runs (Beam->RunType[0] == 'A' or 'a') "
                  "were not enabled at compile time!");
        bail();
#endif
    } else if(rt == 'R') {
        GlobalLog("Internal error, ray run 'R' is not a field mode!");
        bail();
    } else {
        GlobalLog("Invalid Beam->RunType[0] %c!", rt);
        bail();
    }
}

template<bool O3D, bool R3D> inline void RunFieldModesSelInfl(
    bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs, uint32_t cores)
{
    char it = params.Beam->Type[0];
    if(it == 'R') {
#ifdef BHC_INFL_ENABLE_CERVENY_RAYCEN
        if constexpr(!R3D) {
            RunFieldModesSelRun<'R', O3D, R3D>(params, outputs, cores);
        } else {
            GlobalLog("Cerveny ray-centered influence (Beam->Type[0] == 'R') "
                      "is not supported in 3D mode!");
            bail();
        }
#else
        GlobalLog("Cerveny ray-centered influence (Beam->Type[0] == 'R') "
                  "was not enabled at compile time!");
        bail();
#endif
    } else if(it == 'C') {
#ifdef BHC_INFL_ENABLE_CERVENY_CART
        if constexpr(!R3D) {
#ifdef BHC_LIMIT_FEATURES
            if constexpr(!O3D) {
#endif
                RunFieldModesSelRun<'C', O3D, R3D>(params, outputs, cores);
#ifdef BHC_LIMIT_FEATURES
            } else {
                GlobalLog("Nx2D Cerveny Cartesian influence (Beam->Type[0] == 'C') "
                          "is not supported because BHC_LIMIT_FEATURES is enabled!");
                bail();
            }
#endif
        } else {
            GlobalLog("Cerveny Cartesian influence (Beam->Type[0] == 'C') "
                      "is not supported in 3D mode!");
            bail();
        }
#else
        GlobalLog("Cerveny Cartesian influence (Beam->Type[0] == 'C') "
                  "was not enabled at compile time!");
        bail();
#endif
    } else if(it == 'G' || it == '^' || it == ' ') {
#ifdef BHC_INFL_ENABLE_HAT_CART
        RunFieldModesSelRun<'G', O3D, R3D>(params, outputs, cores);
#else
        GlobalLog("Hat Cartesian influence (Beam->Type[0] == 'G', '^', or ' ') "
                  "was not enabled at compile time!");
        bail();
#endif
    } else if(it == 'g') {
#ifdef BHC_INFL_ENABLE_HAT_RAYCEN
        RunFieldModesSelRun<'g', O3D, R3D>(params, outputs, cores);
#else
        GlobalLog("Hat ray-centered influence (Beam->Type[0] == 'g') "
                  "was not enabled at compile time!");
        bail();
#endif
    } else if(it == 'B') {
#ifdef BHC_INFL_ENABLE_GAUSS_CART
        RunFieldModesSelRun<'B', O3D, R3D>(params, outputs, cores);
#else
        GlobalLog("Gaussian Cartesian influence (Beam->Type[0] == 'B') "
                  "was not enabled at compile time!");
        bail();
#endif
    } else if(it == 'b') {
#ifdef BHC_INFL_ENABLE_GAUSS_RAYCEN
#ifdef BHC_LIMIT_FEATURES
        if constexpr(O3D) {
#endif
            RunFieldModesSelRun<'b', O3D, R3D>(params, outputs, cores);
#ifdef BHC_LIMIT_FEATURES
        } else {
            GlobalLog("2D Gaussian RayCen (Beam->Type[0] == 'C') "
                      "is not supported because BHC_LIMIT_FEATURES is enabled!");
            bail();
        }
#endif
#else
        GlobalLog("Gaussian ray-centered influence (Beam->Type[0] == 'b') "
                  "was not enabled at compile time!");
        bail();
#endif
    } else if(it == 'S') {
#ifdef BHC_INFL_ENABLE_SGB
        if constexpr(!R3D) {
            RunFieldModesSelRun<'S', O3D, R3D>(params, outputs, cores);
        } else {
            GlobalLog("Simple Gaussian beams influence (Beam->Type[0] == 'S') "
                      "is not supported in 3D mode!");
        }
#else
        GlobalLog("Simple Gaussian beams influence (Beam->Type[0] == 'S') "
                  "was not enabled at compile time!");
        bail();
#endif
    } else {
        GlobalLog("Invalid Beam->Type[0] %c!", it);
        bail();
    }
}

} // namespace bhc
