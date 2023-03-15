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
#include "common.hpp"
#include "readenv.hpp"
#include "boundary.hpp"
#include "reflect.hpp"
#include "ssp.hpp"
#include "beams.hpp"
#include "run.hpp"

namespace bhc {

#ifdef BHC_BUILD_CUDA
template<bool O3D, bool R3D> void setupGPU(const bhcParams<O3D, R3D> &params)
{
    // Print info about all GPUs and which one is selected
    int num_gpus;
    checkCudaErrors(cudaGetDeviceCount(&num_gpus));
    if(num_gpus <= 0) {
        EXTERR("No CUDA GPUs found; is the driver installed and loaded?");
    }
    int gpuIndex = GetInternal(params)->gpuIndex;
    if(gpuIndex >= num_gpus) {
        EXTERR(
            "You specified CUDA device %d, but there are only %d GPUs", gpuIndex,
            num_gpus);
    }
    cudaDeviceProp cudaProperties;
    for(int g = 0; g < num_gpus; ++g) {
        checkCudaErrors(cudaGetDeviceProperties(&cudaProperties, g));
        if(g == gpuIndex) {
            EXTWARN(
                "CUDA device: %s / compute %d.%d", cudaProperties.name,
                cudaProperties.major, cudaProperties.minor);
        }
        /*
        EXTWARN("%s GPU %d: %s, compute SM %d.%d",
            (g == GetInternal(params)->gpuIndex) ? "-->" : "   "
            g, cudaProperties.name, cudaProperties.major, cudaProperties.minor);
        EXTWARN("      --Global/shared/constant memory: %lli, %d, %d",
            cudaProperties.totalGlobalMem,
            cudaProperties.sharedMemPerBlock,
            cudaProperties.totalConstMem);
        EXTWARN("      --Warp/threads/SMPs: %d, %d, %d" ,
            cudaProperties.warpSize,
            cudaProperties.maxThreadsPerBlock,
            cudaProperties.multiProcessorCount);
        */
    }

    // Store properties about used GPU
    checkCudaErrors(cudaGetDeviceProperties(&cudaProperties, gpuIndex));
    /*
    GetInternal(params)->d_warp       = cudaProperties.warpSize;
    GetInternal(params)->d_maxthreads = cudaProperties.maxThreadsPerBlock;
    */
    GetInternal(params)->d_multiprocs = cudaProperties.multiProcessorCount;
    checkCudaErrors(cudaSetDevice(gpuIndex));
}
#endif

template<bool O3D, bool R3D> bool setup(
    const bhcInit &init, bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs)
{
    try {
        params.internal = new bhcInternal(init);

        Stopwatch sw(GetInternal(params));
        sw.tick();

#ifdef BHC_BUILD_CUDA
        setupGPU(params);
#endif

        // Allocate main structs
        params.Bdry     = nullptr;
        params.bdinfo   = nullptr;
        params.refl     = nullptr;
        params.ssp      = nullptr;
        params.atten    = nullptr;
        params.Pos      = nullptr;
        params.Angles   = nullptr;
        params.freqinfo = nullptr;
        params.Beam     = nullptr;
        params.beaminfo = nullptr;
        outputs.rayinfo = nullptr;
        outputs.eigen   = nullptr;
        outputs.arrinfo = nullptr;
        trackallocate(params, "data structures", params.Bdry);
        trackallocate(params, "data structures", params.bdinfo);
        trackallocate(params, "data structures", params.refl);
        trackallocate(params, "data structures", params.ssp);
        trackallocate(params, "data structures", params.atten);
        trackallocate(params, "data structures", params.Pos);
        trackallocate(params, "data structures", params.Angles);
        trackallocate(params, "data structures", params.freqinfo);
        trackallocate(params, "data structures", params.Beam);
        trackallocate(params, "data structures", params.beaminfo);
        trackallocate(params, "data structures", outputs.rayinfo);
        trackallocate(params, "data structures", outputs.eigen);
        trackallocate(params, "data structures", outputs.arrinfo);

        // Set pointers to null because we always check if they are not null (and
        // deallocate them if so) before allocating them
        // IMPORTANT--if changes are made here, make the same changes in finalize
        params.bdinfo->top.bd          = nullptr;
        params.bdinfo->bot.bd          = nullptr;
        params.refl->bot.r             = nullptr;
        params.refl->top.r             = nullptr;
        params.ssp->cMat               = nullptr;
        params.ssp->czMat              = nullptr;
        params.ssp->Seg.r              = nullptr;
        params.ssp->Seg.x              = nullptr;
        params.ssp->Seg.y              = nullptr;
        params.ssp->Seg.z              = nullptr;
        params.Pos->Sx                 = nullptr;
        params.Pos->Sy                 = nullptr;
        params.Pos->Sz                 = nullptr;
        params.Pos->Rr                 = nullptr;
        params.Pos->Rz                 = nullptr;
        params.Pos->theta              = nullptr;
        params.Pos->t_rcvr             = nullptr;
        params.Angles->alpha.angles    = nullptr;
        params.Angles->beta.angles     = nullptr;
        params.freqinfo->freqVec       = nullptr;
        params.beaminfo->SrcBmPat      = nullptr;
        outputs.rayinfo->results       = nullptr;
        outputs.rayinfo->RayMem        = nullptr;
        outputs.rayinfo->WorkRayMem    = nullptr;
        outputs.uAllSources            = nullptr;
        outputs.eigen->hits            = nullptr;
        outputs.arrinfo->Arr           = nullptr;
        outputs.arrinfo->NArr          = nullptr;
        outputs.arrinfo->MaxNPerSource = nullptr;

        // Fill in default / "constructor" data
        params.fT = RL(1.0e20);
        // Bdry: none
        if constexpr(O3D) {
            params.bdinfo->top.NPts.x = 2;
            params.bdinfo->top.NPts.y = 2;
            params.bdinfo->bot.NPts.x = 2;
            params.bdinfo->bot.NPts.y = 2;
        } else {
            params.bdinfo->top.NPts = 2;
            params.bdinfo->bot.NPts = 2;
        }
        memcpy(params.bdinfo->top.type, "LS", 2);
        memcpy(params.bdinfo->bot.type, "LS", 2);
        // params.refl: none
        params.ssp->dirty      = false;
        params.atten->t        = FL(20.0);
        params.atten->Salinity = FL(35.0);
        params.atten->pH       = FL(8.0);
        params.atten->z_bar    = FL(0.0);
        params.Pos->NSx        = 1;
        params.Pos->NSy        = 1;
        params.Pos->NSz        = 1;
        params.Pos->NRz        = 1;
        params.Pos->NRr        = 1;
        params.Pos->Ntheta     = 1;
        params.Angles->alpha.n = 0;
        params.Angles->beta.n  = 1;
        // LP: not a typo; this is an index, one less than the start of the array,
        // which in Fortran (and in the env file!) is 0. This gets converted to 0-
        // indexed when it is used.
        params.Angles->alpha.iSingle = 0;
        params.Angles->beta.iSingle  = 0;
        params.freqinfo->Nfreq       = 1;
        params.Beam->epsMultiplier   = FL(1.0);
        memcpy(params.Beam->Type, "G S ", 4);
        // params.beaminfo: none
        outputs.rayinfo->RayMemCapacity  = 0;
        outputs.rayinfo->RayMemPoints    = 0;
        outputs.rayinfo->MaxPointsPerRay = 0;
        outputs.rayinfo->NRays           = 0;
        outputs.eigen->neigen            = 0;
        outputs.eigen->memsize           = 0;
        outputs.arrinfo->MaxNArr         = 1;

        ReadEnvironment<O3D, R3D>(params);
        ReadBoundary<O3D, R3D>(
            params, params.Bdry->Top.hs.Opt[4], params.Bdry->Top.hs.Depth,
            &params.bdinfo->top, true); // AlTImetry
        ReadBoundary<O3D, R3D>(
            params, params.Bdry->Bot.hs.Opt[1], params.Bdry->Bot.hs.Depth,
            &params.bdinfo->bot, false);             // BaThYmetry
        ReadReflectionCoefficient<O3D, R3D>(params); // (top and bottom)
        params.beaminfo->SBPFlag = params.Beam->RunType[2];
        ReadPat<O3D, R3D>(params); // Source Beam Pattern
        if constexpr(!O3D) {
            // dummy bearing angles
            params.Pos->Ntheta = 1;
            trackallocate(
                params, "default arrays", params.Pos->theta, params.Pos->Ntheta);
            params.Pos->theta[0] = FL(0.0);
            trackallocate(
                params, "default arrays", params.Pos->t_rcvr, params.Pos->Ntheta);
        }

        // LP: Moved from WriteHeader
        // receiver bearing angles
        if(params.Pos->theta == nullptr) {
            params.Pos->Ntheta = 1;
            trackallocate(params, "default arrays", params.Pos->theta, 1);
            params.Pos->theta[0] = FL(0.0); // dummy bearing angle
            trackallocate(params, "default arrays", params.Pos->t_rcvr, 1);
        }
        // source x-coordinates
        if(params.Pos->Sx == nullptr) {
            trackallocate(params, "default arrays", params.Pos->Sx, 1);
            params.Pos->Sx[0] = FL(0.0); // dummy x-coordinate
            params.Pos->NSx   = 1;
        }
        // source y-coordinates
        if(params.Pos->Sy == nullptr) {
            trackallocate(params, "default arrays", params.Pos->Sy, 1);
            params.Pos->Sy[0] = FL(0.0); // dummy y-coordinate
            params.Pos->NSy   = 1;
        }

        if(params.Beam->deltas == FL(0.0)) {
            // Automatic step size selection
            params.Beam->deltas = (params.Bdry->Bot.hs.Depth - params.Bdry->Top.hs.Depth)
                / FL(10.0);
            if constexpr(!O3D) {
                GetInternal(params)->PRTFile
                    << "\n Step length,       deltas = " << params.Beam->deltas
                    << " m (automatically selected)\n";
            }
        }

        GetInternal(params)->PRTFile << "\n";

        sw.tock("setup");
    } catch(const std::exception &e) {
        EXTWARN("Exception caught in bhc::setup(): %s\n", e.what());
        return false;
    }

    return true;
}

#if BHC_ENABLE_2D
template bool BHC_API setup<false, false>(
    const bhcInit &init, bhcParams<false, false> &params,
    bhcOutputs<false, false> &outputs);
#endif
#if BHC_ENABLE_NX2D
template bool BHC_API setup<true, false>(
    const bhcInit &init, bhcParams<true, false> &params,
    bhcOutputs<true, false> &outputs);
#endif
#if BHC_ENABLE_3D
template bool BHC_API setup<true, true>(
    const bhcInit &init, bhcParams<true, true> &params, bhcOutputs<true, true> &outputs);
#endif

template<bool O3D, bool R3D> void finalize(
    bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs)
{
    // IMPORTANT--if changes are made here, make the same changes in setup
    // (i.e. setting the pointers to nullptr initially)

    trackdeallocate(params, params.bdinfo->top.bd);
    trackdeallocate(params, params.bdinfo->bot.bd);
    trackdeallocate(params, params.refl->bot.r);
    trackdeallocate(params, params.refl->top.r);
    trackdeallocate(params, params.ssp->cMat);
    trackdeallocate(params, params.ssp->czMat);
    trackdeallocate(params, params.ssp->Seg.r);
    trackdeallocate(params, params.ssp->Seg.x);
    trackdeallocate(params, params.ssp->Seg.y);
    trackdeallocate(params, params.ssp->Seg.z);
    trackdeallocate(params, params.Pos->iSz);
    trackdeallocate(params, params.Pos->iRz);
    trackdeallocate(params, params.Pos->Sx);
    trackdeallocate(params, params.Pos->Sy);
    trackdeallocate(params, params.Pos->Sz);
    trackdeallocate(params, params.Pos->Rr);
    trackdeallocate(params, params.Pos->Rz);
    trackdeallocate(params, params.Pos->ws);
    trackdeallocate(params, params.Pos->wr);
    trackdeallocate(params, params.Pos->theta);
    trackdeallocate(params, params.Pos->t_rcvr);
    trackdeallocate(params, params.Angles->alpha.angles);
    trackdeallocate(params, params.Angles->beta.angles);
    trackdeallocate(params, params.freqinfo->freqVec);
    trackdeallocate(params, params.beaminfo->SrcBmPat);
    trackdeallocate(params, outputs.rayinfo->results);
    trackdeallocate(params, outputs.rayinfo->RayMem);
    trackdeallocate(params, outputs.rayinfo->WorkRayMem);
    trackdeallocate(params, outputs.uAllSources);
    trackdeallocate(params, outputs.eigen->hits);
    trackdeallocate(params, outputs.arrinfo->Arr);
    trackdeallocate(params, outputs.arrinfo->NArr);
    trackdeallocate(params, outputs.arrinfo->MaxNPerSource);

    trackdeallocate(params, params.Bdry);
    trackdeallocate(params, params.bdinfo);
    trackdeallocate(params, params.refl);
    trackdeallocate(params, params.ssp);
    trackdeallocate(params, params.atten);
    trackdeallocate(params, params.Pos);
    trackdeallocate(params, params.Angles);
    trackdeallocate(params, params.freqinfo);
    trackdeallocate(params, params.Beam);
    trackdeallocate(params, params.beaminfo);
    trackdeallocate(params, outputs.rayinfo);
    trackdeallocate(params, outputs.eigen);
    trackdeallocate(params, outputs.arrinfo);

    if(GetInternal(params)->usedMemory != 0) {
        EXTWARN(
            "Amount of memory leaked: %" PRIu64 " bytes",
            GetInternal(params)->usedMemory);
    }

    delete GetInternal(params);
    params.internal = nullptr;
}

#if BHC_ENABLE_2D
template void BHC_API finalize<false, false>(
    bhcParams<false, false> &params, bhcOutputs<false, false> &outputs);
#endif
#if BHC_ENABLE_NX2D
template void BHC_API
finalize<true, false>(bhcParams<true, false> &params, bhcOutputs<true, false> &outputs);
#endif
#if BHC_ENABLE_3D
template void BHC_API
finalize<true, true>(bhcParams<true, true> &params, bhcOutputs<true, true> &outputs);
#endif

} // namespace bhc
