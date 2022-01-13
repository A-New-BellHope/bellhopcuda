#include "setup.hpp"
#include "main.hpp"
#include "output.hpp"

std::ofstream PRTFile, ARRFile;
LDOFile RAYFile;
DirectOFile SHDFile;
std::string Title;
real fT;
BdryType *Bdry;
BdryInfo *bdinfo;
ReflectionInfo *refl;
SSPStructure *ssp;
AttenInfo *atten;
Position *Pos;
AnglesStructure *Angles;
FreqInfo *freqinfo;
BeamStructure *Beam;
BeamInfo *beaminfo;

cpxf *uAllSources;

__global__ void 
__launch_bounds__(512, 1)
TLModeKernel(cpxf *uAllSources_, 
    const BdryType *ConstBdry_, const BdryInfo *bdinfo_, const ReflectionInfo *refl_,
    const SSPStructure *ssp_, const Position *Pos_, const AnglesStructure *Angles_,
    const FreqInfo *freqinfo_, const BeamStructure *Beam_, const BeamInfo *beaminfo_)
{
    for(int32_t job = blockIdx.x * blockDim.x + threadIdx.x; ; job += gridDim.x * blockDim.x){
        int32_t isrc, ialpha;
        if(!GetJobIndices(isrc, ialpha, job, Pos_, Angles_)) break;
        
        real SrcDeclAngle;
        MainTLMode(isrc, ialpha, SrcDeclAngle, uAllSources_,
            ConstBdry_, bdinfo_, refl_, ssp_, Pos_, Angles_, freqinfo_, Beam_, beaminfo_);
    }
}

int m_gpu, d_warp, d_maxthreads, d_multiprocs;
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

int main(int argc, char **argv)
{
    if(argc != 2){
        std::cout << "Must provide FileRoot as command-line parameter\n";
        std::abort();
    }
    std::string FileRoot = argv[1];
    
    setupGPU();
    setup(FileRoot, PRTFile, RAYFile, ARRFile, SHDFile, Title, fT,
        Bdry, bdinfo, refl, ssp, atten, Pos, Angles, freqinfo, Beam, beaminfo);   
    
    if(Beam->RunType[0] == 'R'){
        std::cout << "Ray runs not implemented in CUDA\n";
        std::abort();
    }else if(Beam->RunType[0] == 'C' || Beam->RunType[0] == 'S' || Beam->RunType[0] == 'I'){
        // TL mode
        InitTLMode(uAllSources, Pos, Beam);
        
        Stopwatch sw;
        sw.tick();
        TLModeKernel<<<d_multiprocs,512>>>(uAllSources, 
            Bdry, bdinfo, refl, ssp, Pos, Angles, freqinfo, Beam, beaminfo);
        syncAndCheckKernelErrors("TLModeKernel");
        sw.tock();
        
        //std::cout << "Output\n";
        FinalizeTLMode(uAllSources, SHDFile, ssp, Pos, Angles, freqinfo, Beam);
    }else{
        std::cout << "Not yet implemented RunType " << Beam->RunType[0] << "\n";
        std::abort();
    }
}
