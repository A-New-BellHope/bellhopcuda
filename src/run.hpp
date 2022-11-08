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
extern int m_gpu = 0, d_warp, d_maxthreads, d_multiprocs;
#endif

template<bool O3D, bool R3D> inline void InitSelectedMode(
    bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs, bool singlethread)
{
    // this is always called from run_* so update intermediate params
    PrintFileEmu &PRTFile = *(PrintFileEmu *)params.internal;
    bhc::UpdateSSP(params.freqinfo->freq0, params.fT, params.ssp, PRTFile, params.atten);

    // Common
    int32_t ns = params.Pos->NSx * params.Pos->NSy * params.Pos->NSz;
    BASSERT(ns >= 1);

    // irregular or rectilinear grid
    params.Pos->NRz_per_range = IsIrregularGrid(params.Beam) ? 1 : params.Pos->NRz;

    // Mode specific
    if(IsRayRun(params.Beam)) {
        InitRayMode<O3D, R3D>(outputs.rayinfo, params);
    } else if(IsTLRun(params.Beam)) {
        InitTLMode(outputs.uAllSources, params.Pos);
    } else if(IsEigenraysRun(params.Beam)) {
        InitEigenMode(outputs.eigen);
    } else if(IsArrivalsRun(params.Beam)) {
        InitArrivalsMode(
            outputs.arrinfo, singlethread, R3D, params.Pos,
            *(PrintFileEmu *)params.internal);
    } else {
        GlobalLog("Invalid RunType %c\n", params.Beam->RunType[0]);
        std::abort();
    }
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
        RunFieldModesImpl<CfgSel<RT, IT, 'N'>, O3D, R3D>(params, outputs, cores);
    } else if(st == 'C') {
        RunFieldModesImpl<CfgSel<RT, IT, 'C'>, O3D, R3D>(params, outputs, cores);
    } else if(st == 'S') {
        RunFieldModesImpl<CfgSel<RT, IT, 'S'>, O3D, R3D>(params, outputs, cores);
    } else if(st == 'P') {
        RunFieldModesImpl<CfgSel<RT, IT, 'P'>, O3D, R3D>(params, outputs, cores);
    } else if(st == 'Q') {
        RunFieldModesImpl<CfgSel<RT, IT, 'Q'>, O3D, R3D>(params, outputs, cores);
    } else if(st == 'H') {
        RunFieldModesImpl<CfgSel<RT, IT, 'H'>, O3D, R3D>(params, outputs, cores);
    } else if(st == 'A') {
        RunFieldModesImpl<CfgSel<RT, IT, 'A'>, O3D, R3D>(params, outputs, cores);
    } else {
        GlobalLog("Invalid ssp->Type %c for field modes!", st);
        bail();
    }
}

template<char IT, bool O3D, bool R3D> inline void RunFieldModesSelRun(
    bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs, uint32_t cores)
{
    char rt = params.Beam->RunType[0];
    if(rt == 'C' || rt == 'S' || rt == 'I') {
        RunFieldModesSelSSP<'C', IT, O3D, R3D>(params, outputs, cores);
    } else if(rt == 'E') {
        RunFieldModesSelSSP<'E', IT, O3D, R3D>(params, outputs, cores);
    } else if(rt == 'A' || rt == 'a') {
        RunFieldModesSelSSP<'A', IT, O3D, R3D>(params, outputs, cores);
    } else {
        GlobalLog("Invalid RunType %c for field modes!", rt);
        bail();
    }
}

template<bool O3D, bool R3D> inline void RunFieldModesSelInfl(
    bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs, uint32_t cores)
{
    char it = params.Beam->Type[0];
    if(it == 'R') {
        RunFieldModesSelRun<'R', O3D, R3D>(params, outputs, cores);
    } else if(it == 'C') {
        RunFieldModesSelRun<'C', O3D, R3D>(params, outputs, cores);
    } else if(it == 'G' || it == '^' || it == ' ') {
        RunFieldModesSelRun<'G', O3D, R3D>(params, outputs, cores);
    } else if(it == 'g') {
        RunFieldModesSelRun<'g', O3D, R3D>(params, outputs, cores);
    } else if(it == 'B') {
        RunFieldModesSelRun<'B', O3D, R3D>(params, outputs, cores);
    } else if(it == 'b') {
        RunFieldModesSelRun<'b', O3D, R3D>(params, outputs, cores);
    } else if(it == 'S') {
        RunFieldModesSelRun<'S', O3D, R3D>(params, outputs, cores);
    } else {
        GlobalLog("Invalid Beam type %c for field modes!", it);
        bail();
    }
}

} // namespace bhc
