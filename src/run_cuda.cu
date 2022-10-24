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
#include "run.hpp"

namespace bhc {

#define NUM_THREADS 256

template<bool O3D, bool R3D> __global__ void __launch_bounds__(NUM_THREADS, 1)
FieldModesKernel(bhcParams<O3D, R3D> params, bhcOutputs<O3D, R3D> outputs)
{
    for(int32_t job = blockIdx.x * blockDim.x + threadIdx.x; ; job += gridDim.x * blockDim.x){
        RayInitInfo rinit;
        if(!GetJobIndices<O3D>(rinit, job, params.Pos, params.Angles)) break;
        
        MainFieldModes<O3D, R3D>(rinit, outputs.uAllSources,
            params.Bdry, params.bdinfo, params.refl, params.ssp, params.Pos,
            params.Angles, params.freqinfo, params.Beam, params.beaminfo, 
            outputs.eigen, outputs.arrinfo);
    }
}

#if BHC_ENABLE_2D
template __global__ void __launch_bounds__(NUM_THREADS, 1) FieldModesKernel<false, false>(
    bhcParams<false, false> params, bhcOutputs<false, false> outputs);
#endif
#if BHC_ENABLE_NX2D
template __global__ void __launch_bounds__(NUM_THREADS, 1) FieldModesKernel<true, false>(
    bhcParams<true, false> params, bhcOutputs<true, false> outputs);
#endif
#if BHC_ENABLE_3D
template __global__ void __launch_bounds__(NUM_THREADS, 1) FieldModesKernel<true, true>(
    bhcParams<true, true> params, bhcOutputs<true, true> outputs);
#endif

int m_gpu = 0, d_warp, d_maxthreads, d_multiprocs;
void setupGPU()
{
    //Print info about all GPUs and which one is selected
    int num_gpus;
    checkCudaErrors(cudaGetDeviceCount(&num_gpus));
    BASSERT(num_gpus >= 1);
    cudaDeviceProp cudaProperties;
    for(int g=0; g<num_gpus; ++g){
        checkCudaErrors(cudaGetDeviceProperties(&cudaProperties, g));
        if(g == m_gpu){
            GlobalLog("CUDA device: %s / compute %d.%d\n",
                cudaProperties.name, cudaProperties.major, cudaProperties.minor);
        }
        /*
        GlobalLog("%s", (g == m_gpu) ? "--> " : "    ");
        GlobalLog("GPU %d: %s, compute SM %d.%d\n",
            g, cudaProperties.name, cudaProperties.major, cudaProperties.minor);
        GlobalLog("      --Global/shared/constant memory: %lli, %d, %d\n",
            cudaProperties.totalGlobalMem,
            cudaProperties.sharedMemPerBlock,
            cudaProperties.totalConstMem);
        GlobalLog("      --Warp/threads/SMPs: %d, %d, %d\n" ,
            cudaProperties.warpSize,
            cudaProperties.maxThreadsPerBlock,
            cudaProperties.multiProcessorCount);
        */
    }
    
    //Store properties about used GPU
    checkCudaErrors(cudaGetDeviceProperties(&cudaProperties, m_gpu));
    d_warp = cudaProperties.warpSize;
    d_maxthreads = cudaProperties.maxThreadsPerBlock;
    d_multiprocs = cudaProperties.multiProcessorCount;
    checkCudaErrors(cudaSetDevice(m_gpu));
}

template<bool O3D, bool R3D> bool run_cuda(
    bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs)
{
    if(!api_okay) return false;
    
    try{
    
    InitSelectedMode<O3D, R3D>(params, outputs, false);
    FieldModesKernel<O3D, R3D><<<d_multiprocs,NUM_THREADS>>>(params, outputs);
    syncAndCheckKernelErrors("FieldModesKernel");
    
    }catch(const std::exception &e){
        api_okay = false;
        PrintFileEmu &PRTFile = *(PrintFileEmu*)params.internal;
        PRTFile << "Exception caught:\n" << e.what() << "\n";
    }
    
    return api_okay;
}

#if BHC_ENABLE_2D
template bool run_cuda<false, false>(
    bhcParams<false, false> &params, bhcOutputs<false, false> &outputs);
#endif
#if BHC_ENABLE_NX2D
template bool run_cuda<true, false>(
    bhcParams<true, false> &params, bhcOutputs<true, false> &outputs);
#endif
#if BHC_ENABLE_3D
template bool run_cuda<true, true>(
    bhcParams<true, true> &params, bhcOutputs<true, true> &outputs);
#endif

template<bool O3D, bool R3D> bool run(
    bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs, bool singlethread)
{
    if(singlethread){
        GlobalLog("Single threaded mode is nonsense on CUDA, ignoring\n");
    }
    if(IsRayRun(params.Beam)){
        return run_cxx<O3D, R3D>(params, outputs, false);
    }else{
        return run_cuda<O3D, R3D>(params, outputs);
    }
}

#if BHC_ENABLE_2D
template bool BHC_API run<false, false>(
    bhcParams<false, false> &params, bhcOutputs<false, false> &outputs, bool singlethread);
#endif
#if BHC_ENABLE_NX2D
template bool BHC_API run<true, false>(
    bhcParams<true, false> &params, bhcOutputs<true, false> &outputs, bool singlethread);
#endif
#if BHC_ENABLE_3D
template bool BHC_API run<true, true>(
    bhcParams<true, true> &params, bhcOutputs<true, true> &outputs, bool singlethread); 
#endif


}
