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

#define MAX_THREADS 256
__global__ void 
__launch_bounds__(MAX_THREADS, 1)
FieldModesKernel(bhcParams params, bhcOutputs outputs)
{
    for(int32_t job = blockIdx.x * blockDim.x + threadIdx.x; ; job += gridDim.x * blockDim.x){
        int32_t isrc, ialpha;
        if(!GetJobIndices(isrc, ialpha, job, params.Pos, params.Angles)) break;
        
        real SrcDeclAngle;
        MainFieldModes(isrc, ialpha, SrcDeclAngle, outputs.uAllSources,
            params.Bdry, params.bdinfo, params.refl, params.ssp, params.Pos,
            params.Angles, params.freqinfo, params.Beam, params.beaminfo, 
            outputs.eigen, outputs.arrinfo);
    }
}

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
            std::cout << "CUDA device: " << cudaProperties.name << " / compute "
                << cudaProperties.major << "." << cudaProperties.minor << "\n";
        }
        /*
        std::cout << ((g == m_gpu) ? "--> " : "    ");
        std::cout << "GPU " << g << ": " << cudaProperties.name << ", compute SM " 
            << cudaProperties.major << "." << cudaProperties.minor << "\n";
        std::cout << "      --Global/shared/constant memory: " 
            << cudaProperties.totalGlobalMem << ", " 
            << cudaProperties.sharedMemPerBlock << ", "
            << cudaProperties.totalConstMem << "\n";
        std::cout << "      --Warp/threads/SMPs: " 
            << cudaProperties.warpSize << ", "
            << cudaProperties.maxThreadsPerBlock << ", "
            << cudaProperties.multiProcessorCount << "\n";
        */
    }
    
    //Store properties about used GPU
    checkCudaErrors(cudaGetDeviceProperties(&cudaProperties, m_gpu));
    d_warp = cudaProperties.warpSize;
    d_maxthreads = cudaProperties.maxThreadsPerBlock;
    d_multiprocs = cudaProperties.multiProcessorCount;
    checkCudaErrors(cudaSetDevice(m_gpu));
}

bool run_cuda(const bhcParams &params, bhcOutputs &outputs)
{
    if(!api_okay) return false;
    
    try{
    
    InitSelectedMode(params, outputs, false);
    FieldModesKernel<<<d_multiprocs,MAX_THREADS>>>(params, outputs);
    syncAndCheckKernelErrors("FieldModesKernel");
    
    }catch(const std::exception &e){
        api_okay = false;
        PrintFileEmu &PRTFile = *(PrintFileEmu*)params.internal;
        PRTFile << "Exception caught:\n" << e.what() << "\n";
    }
    
    return api_okay;
}

BHC_API bool run(const bhcParams &params, bhcOutputs &outputs, bool singlethread)
{
    if(singlethread){
        std::cout << "Single threaded mode is nonsense on CUDA, ignoring\n";
    }
    if(params.Beam->RunType[0] == 'R'){
        return run_cxx(params, outputs, false);
    }else{
        return run_cuda(params, outputs);
    }
}

}
