#include "run.hpp"

__global__ void 
__launch_bounds__(512, 1)
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

void run_cuda(std::ostream &PRTFile, const bhcParams &params, bhcOutputs &outputs)
{
    InitSelectedMode(PRTFile, params, outputs, false);
    FieldModesKernel<<<d_multiprocs,512>>>(params, outputs);
    syncAndCheckKernelErrors("FieldModesKernel");
}

BHC_API void run(std::ostream &PRTFile, const bhcParams &params, bhcOutputs &outputs,
    bool singlethread)
{
    if(singlethread){
        std::cout << "Single threaded mode is nonsense on CUDA, ignoring\n";
    }
    if(params.Beam->RunType[0] == 'R'){
        run_cxx(PRTFile, params, outputs, false);
    }else{
        run_cuda(PRTFile, params, outputs);
    }
}
