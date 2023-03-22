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

        outputs.rayinfo->results       = nullptr;
        outputs.rayinfo->RayMem        = nullptr;
        outputs.rayinfo->WorkRayMem    = nullptr;
        outputs.uAllSources            = nullptr;
        outputs.eigen->hits            = nullptr;
        outputs.arrinfo->Arr           = nullptr;
        outputs.arrinfo->NArr          = nullptr;
        outputs.arrinfo->MaxNPerSource = nullptr;

        // Fill in default / "constructor" data
        // Bdry: none
        // params.refl: none

        // params.beaminfo: none
        outputs.rayinfo->RayMemCapacity  = 0;
        outputs.rayinfo->RayMemPoints    = 0;
        outputs.rayinfo->MaxPointsPerRay = 0;
        outputs.rayinfo->NRays           = 0;
        outputs.eigen->neigen            = 0;
        outputs.eigen->memsize           = 0;
        outputs.arrinfo->MaxNArr         = 1;

        ReadEnvironment<O3D, R3D>(params);

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
