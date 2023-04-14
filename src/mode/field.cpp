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
#include "field.hpp"
#include "../common_run.hpp"

namespace bhc { namespace mode {

template<char RT, char IT, bool O3D, bool R3D> inline void RunFieldModesSelSSP(
    bhcParams<O3D> &params, bhcOutputs<O3D, R3D> &outputs)
{
    char st = params.ssp->Type;
    if(st == 'N') {
#ifdef BHC_SSP_ENABLE_N2LINEAR
        RunFieldModesImpl<CfgSel<RT, IT, 'N'>, O3D, R3D>(params, outputs);
#else
        EXTERR("N2-linear SSP (ssp->Type == 'N') was not enabled at compile time!");
#endif
    } else if(st == 'C') {
#ifdef BHC_SSP_ENABLE_CLINEAR
        RunFieldModesImpl<CfgSel<RT, IT, 'C'>, O3D, R3D>(params, outputs);
#else
        EXTERR("C-linear SSP (ssp->Type == 'C') was not enabled at compile time!");
#endif
    } else if(st == 'S') {
#ifdef BHC_SSP_ENABLE_CUBIC
        RunFieldModesImpl<CfgSel<RT, IT, 'S'>, O3D, R3D>(params, outputs);
#else
        EXTERR("Cubic spline SSP (ssp->Type == 'S') was not enabled at compile time!");
#endif
    } else if(st == 'P') {
#ifdef BHC_SSP_ENABLE_PCHIP
#ifdef BHC_LIMIT_FEATURES
        if constexpr(!O3D) {
#endif
            RunFieldModesImpl<CfgSel<RT, IT, 'P'>, O3D, R3D>(params, outputs);
#ifdef BHC_LIMIT_FEATURES
        } else {
            EXTERR("Nx2D or 3D PCHIP SSP not supported"
                   "because BHC_LIMIT_FEATURES enabled!");
        }
#endif
#else
        EXTERR("PCHIP SSP (ssp->Type == 'P') was not enabled at compile time!");
#endif
    } else if(st == 'Q') {
#ifdef BHC_SSP_ENABLE_QUAD
        if constexpr(!O3D) {
            RunFieldModesImpl<CfgSel<RT, IT, 'Q'>, O3D, R3D>(params, outputs);
        } else {
            EXTERR("Quad SSP not supported in Nx2D or 3D mode!");
        }
#else
        EXTERR("Quad SSP (ssp->Type == 'Q') was not enabled at compile time!");
#endif
    } else if(st == 'H') {
#ifdef BHC_SSP_ENABLE_HEXAHEDRAL
        if constexpr(O3D) {
            RunFieldModesImpl<CfgSel<RT, IT, 'H'>, O3D, R3D>(params, outputs);
        } else {
            EXTERR("Hexahedral SSP not supported in 2D mode!");
        }
#else
        EXTERR("Hexahedral SSP (ssp->Type == 'H') was not enabled at compile time!");
#endif
    } else if(st == 'A') {
#ifdef BHC_SSP_ENABLE_ANALYTIC
        RunFieldModesImpl<CfgSel<RT, IT, 'A'>, O3D, R3D>(params, outputs);
#else
        EXTERR("Analytic SSP (ssp->Type == 'A') was not enabled at compile time!");
#endif
    } else {
        EXTERR("Invalid ssp->Type %c!", st);
    }
}

template<char IT, bool O3D, bool R3D> inline void RunFieldModesSelRun(
    bhcParams<O3D> &params, bhcOutputs<O3D, R3D> &outputs)
{
    char rt = params.Beam->RunType[0];
    if(rt == 'C' || rt == 'S' || rt == 'I') {
#ifdef BHC_RUN_ENABLE_TL
        RunFieldModesSelSSP<'C', IT, O3D, R3D>(params, outputs);
#else
        EXTERR("Transmission loss runs (Beam->RunType[0] == 'C', 'S', or 'I') "
               "were not enabled at compile time!");
#endif
    } else if(rt == 'E') {
#ifdef BHC_RUN_ENABLE_EIGENRAYS
        if constexpr(InflType<IT>::IsCerveny()) {
            EXTERR("Cerveny influence does not support eigenrays!");
        } else {
            RunFieldModesSelSSP<'E', IT, O3D, R3D>(params, outputs);
        }
#else
        EXTERR("Eigenrays runs (Beam->RunType[0] == 'E') "
               "were not enabled at compile time!");
#endif
    } else if(rt == 'A' || rt == 'a') {
#ifdef BHC_RUN_ENABLE_ARRIVALS
        if constexpr(InflType<IT>::IsCerveny()) {
            EXTERR("Cerveny influence does not support arrivals!");
        } else {
            RunFieldModesSelSSP<'A', IT, O3D, R3D>(params, outputs);
        }
#else
        EXTERR("Arrivals runs (Beam->RunType[0] == 'A' or 'a') "
               "were not enabled at compile time!");
#endif
    } else if(rt == 'R') {
        EXTERR("Internal error, ray run 'R' is not a field mode!");
    } else {
        EXTERR("Invalid Beam->RunType[0] %c!", rt);
    }
}

template<bool O3D, bool R3D> inline void RunFieldModesSelInfl(
    bhcParams<O3D> &params, bhcOutputs<O3D, R3D> &outputs)
{
#ifdef BHC_BUILD_CUDA
    checkCudaErrors(cudaSetDevice(GetInternal(params)->gpuIndex));
#endif
    char it = params.Beam->Type[0];
    if(it == 'R') {
#ifdef BHC_INFL_ENABLE_CERVENY_RAYCEN
        if constexpr(!R3D) {
            RunFieldModesSelRun<'R', O3D, R3D>(params, outputs);
        } else {
            EXTERR("Cerveny ray-centered influence (Beam->Type[0] == 'R') "
                   "is not supported in 3D mode!");
        }
#else
        EXTERR("Cerveny ray-centered influence (Beam->Type[0] == 'R') "
               "was not enabled at compile time!");
#endif
    } else if(it == 'C') {
#ifdef BHC_INFL_ENABLE_CERVENY_CART
        if constexpr(!R3D) {
#ifdef BHC_LIMIT_FEATURES
            if constexpr(!O3D) {
#endif
                RunFieldModesSelRun<'C', O3D, R3D>(params, outputs);
#ifdef BHC_LIMIT_FEATURES
            } else {
                EXTERR("Nx2D Cerveny Cartesian influence (Beam->Type[0] == 'C') "
                       "is not supported because BHC_LIMIT_FEATURES is enabled!");
            }
#endif
        } else {
            EXTERR("Cerveny Cartesian influence (Beam->Type[0] == 'C') "
                   "is not supported in 3D mode!");
        }
#else
        EXTERR("Cerveny Cartesian influence (Beam->Type[0] == 'C') "
               "was not enabled at compile time!");
#endif
    } else if(it == 'G' || it == '^' || it == ' ' || it == 'B') {
#ifdef BHC_INFL_ENABLE_GEOM_CART
        RunFieldModesSelRun<'G', O3D, R3D>(params, outputs);
#else
        EXTERR("Geometric Cartesian influence (Beam->Type[0] == 'G', '^', ' ' "
               "hat / 'B' Gaussian) was not enabled at compile time!");
#endif
    } else if(it == 'g' || it == 'b') {
#ifdef BHC_INFL_ENABLE_GEOM_RAYCEN
#ifdef BHC_LIMIT_FEATURES
        if(it == 'b') {
            if constexpr(!O3D) {
                EXTERR("2D Gaussian RayCen (Beam->Type[0] == 'b') "
                       "is not supported because BHC_LIMIT_FEATURES is enabled!");
            }
        }
#endif
        RunFieldModesSelRun<'g', O3D, R3D>(params, outputs);
#else
        EXTERR("Geometric ray-centered influence (Beam->Type[0] == 'g' hat / "
               "'b' Gaussian) was not enabled at compile time!");
#endif
    } else if(it == 'S') {
#ifdef BHC_INFL_ENABLE_SGB
        if constexpr(!R3D) {
            RunFieldModesSelRun<'S', O3D, R3D>(params, outputs);
        } else {
            EXTERR("Simple Gaussian beams influence (Beam->Type[0] == 'S') "
                   "is not supported in 3D mode!");
        }
#else
        EXTERR("Simple Gaussian beams influence (Beam->Type[0] == 'S') "
               "was not enabled at compile time!");
#endif
    } else {
        EXTERR("Invalid Beam->Type[0] %c!", it);
    }
}

#if BHC_ENABLE_2D
template void RunFieldModesSelInfl<false, false>(
    bhcParams<false> &params, bhcOutputs<false, false> &outputs);
#endif
#if BHC_ENABLE_NX2D
template void RunFieldModesSelInfl<true, false>(
    bhcParams<true> &params, bhcOutputs<true, false> &outputs);
#endif
#if BHC_ENABLE_3D
template void RunFieldModesSelInfl<true, true>(
    bhcParams<true> &params, bhcOutputs<true, true> &outputs);
#endif

}} // namespace bhc::mode
