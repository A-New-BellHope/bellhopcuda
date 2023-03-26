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
#include "common_setup.hpp"

#include "module/paramsmodule.hpp"
#include "module/atten.hpp"
#include "module/title.hpp"
#include "module/freq0.hpp"
#include "module/nmedia.hpp"
#include "module/topopt.hpp"
#include "module/boundarycond.hpp"
#include "module/ssp.hpp"
#include "module/botopt.hpp"
#include "module/sxsy.hpp"
#include "module/szrz.hpp"
#include "module/rcvrranges.hpp"
#include "module/rcvrbearings.hpp"
#include "module/freqvec.hpp"
#include "module/runtype.hpp"
#include "module/rayangles.hpp"
#include "module/beaminfo.hpp"
#include "module/boundary.hpp"
#include "module/reflcoef.hpp"
#include "module/pat.hpp"

#include "mode/modemodule.hpp"
#include "mode/ray.hpp"
#include "mode/field.hpp"
#include "mode/tl.hpp"
#include "mode/eigen.hpp"
#include "mode/arr.hpp"

namespace bhc {

namespace module {

template<bool O3D, bool R3D> class ModulesList {
public:
    ModulesList()
    {
        modules.push_back(new Atten<O3D, R3D>());
        modules.push_back(new Title<O3D, R3D>());
        modules.push_back(new Freq0<O3D, R3D>());
        modules.push_back(new NMedia<O3D, R3D>());
        modules.push_back(new TopOpt<O3D, R3D>());
        modules.push_back(new BoundaryCondTop<O3D, R3D>());
        modules.push_back(new SSP<O3D, R3D>());
        modules.push_back(new BotOpt<O3D, R3D>());
        modules.push_back(new BoundaryCondBot<O3D, R3D>());
        modules.push_back(new SxSy<O3D, R3D>());
        modules.push_back(new SzRz<O3D, R3D>());
        modules.push_back(new RcvrRanges<O3D, R3D>());
        modules.push_back(new RcvrBearings<O3D, R3D>());
        modules.push_back(new FreqVec<O3D, R3D>());
        modules.push_back(new RunType<O3D, R3D>());
        modules.push_back(new RayAnglesElevation<O3D, R3D>());
        modules.push_back(new RayAnglesBearing<O3D, R3D>());
        modules.push_back(new BeamInfo<O3D, R3D>());
        modules.push_back(new Altimetry<O3D, R3D>());
        modules.push_back(new Bathymetry<O3D, R3D>());
        modules.push_back(new BRC<O3D, R3D>());
        modules.push_back(new TRC<O3D, R3D>());
        modules.push_back(new Pat<O3D, R3D>());
    }
    ~ModulesList()
    {
        for(auto *module : modules) delete module;
    }

    const std::vector<ParamsModule<O3D, R3D> *> &list() const { return modules; }

private:
    std::vector<ParamsModule<O3D, R3D> *> modules;
};

#if BHC_ENABLE_2D
template class ParamsModule<false, false>;
template class ModulesList<false, false>;
#endif
#if BHC_ENABLE_NX2D
template class ParamsModule<true, false>;
template class ModulesList<true, false>;
#endif
#if BHC_ENABLE_3D
template class ParamsModule<true, true>;
template class ModulesList<true, true>;
#endif

} // namespace module

namespace mode {

template<bool O3D, bool R3D> class ModesList {
public:
    ModesList()
    {
        modes.push_back(new Ray<O3D, R3D>());
        modes.push_back(new TL<O3D, R3D>());
        modes.push_back(new Eigen<O3D, R3D>());
        modes.push_back(new Arr<O3D, R3D>());
    }
    ~ModesList()
    {
        for(auto *mode : modes) delete mode;
    }

    const std::vector<ModeModule<O3D, R3D> *> &list() const { return modes; }

private:
    std::vector<ModeModule<O3D, R3D> *> modes;
};

#if BHC_ENABLE_2D
template class ModeModule<false, false>;
template class Field<false, false>;
template class ModesList<false, false>;
#endif
#if BHC_ENABLE_NX2D
template class ModeModule<true, false>;
template class Field<true, false>;
template class ModesList<true, false>;
#endif
#if BHC_ENABLE_3D
template class ModeModule<true, true>;
template class Field<true, true>;
template class ModesList<true, true>;
#endif

} // namespace mode

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

////////////////////////////////////////////////////////////////////////////////

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
        params.sbp      = nullptr;
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
        trackallocate(params, "data structures", params.sbp);
        trackallocate(params, "data structures", outputs.rayinfo);
        trackallocate(params, "data structures", outputs.eigen);
        trackallocate(params, "data structures", outputs.arrinfo);

        module::ModulesList<O3D, R3D> modules;
        mode::ModesList<O3D, R3D> modes;
        for(auto *m : modules.list()) m->Init(params);
        for(auto *m : modes.list()) m->Init(outputs);

        if(GetInternal(params)->noEnvFil) {
            for(auto *m : modules.list()) {
                m->SetupPre(params);
                m->Default(params);
                m->SetupPost(params);
            }
        } else {
            // Values only initialized once--reused from top to ssp, and ssp to bot
            HSInfo RecycledHS;
            RecycledHS.alphaR = FL(1500.0);
            RecycledHS.betaR  = FL(0.0);
            RecycledHS.rho    = FL(1.0);
            RecycledHS.alphaI = FL(0.0);
            RecycledHS.betaI  = FL(0.0);

            PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;
            PRTFile << BHC_PROGRAMNAME << (R3D ? "3D" : O3D ? "Nx2D" : "") << "\n\n";

            // Open the environmental file
            LDIFile ENVFile(GetInternal(params), GetInternal(params)->FileRoot + ".env");
            if(!ENVFile.Good()) {
                PRTFile << "ENVFile = " << GetInternal(params)->FileRoot << ".env\n";
                EXTERR(BHC_PROGRAMNAME
                       " - ReadEnvironment: Unable to open the environmental file");
            }

            for(auto *m : modules.list()) {
                m->SetupPre(params);
                m->Read(params, ENVFile, RecycledHS);
                m->SetupPost(params);
            }
        }
        for(auto *m : modules.list()) {
            m->Validate(params);
            m->Echo(params);
        }

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

////////////////////////////////////////////////////////////////////////////////

template<bool O3D, bool R3D> inline mode::ModeModule<O3D, R3D> *GetMode(
    const bhcParams<O3D, R3D> &params)
{
    if(IsRayRun(params.Beam)) {
        return new mode::Ray<O3D, R3D>();
    } else if(IsTLRun(params.Beam)) {
        return new mode::TL<O3D, R3D>();
    } else if(IsEigenraysRun(params.Beam)) {
        return new mode::Eigen<O3D, R3D>();
    } else if(IsArrivalsRun(params.Beam)) {
        return new mode::Arr<O3D, R3D>();
    } else {
        EXTERR("Invalid RunType %c\n", params.Beam->RunType[0]);
        return nullptr;
    }
}

template<bool O3D, bool R3D> bool run(
    bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs)
{
    try {
        Stopwatch sw(GetInternal(params));

        sw.tick();
        module::ModulesList<O3D, R3D> modules;
        for(auto *m : modules.list()) m->Validate(params);
        for(auto *m : modules.list()) m->Preprocess(params);
        auto *mo = GetMode<O3D, R3D>(params);
        mo->Preprocess(params, outputs);
        sw.tock("Preprocess");

        sw.tick();
        mo->Run(params, outputs);
        sw.tock("Run");

        sw.tick();
        mo->Postprocess(params, outputs);
        sw.tock("Postprocess");

        delete mo;
    } catch(const std::exception &e) {
        EXTWARN("Exception caught in bhc::run(): %s\n", e.what());
        return false;
    }

    return true;
}

#if BHC_ENABLE_2D
template bool BHC_API
run<false, false>(bhcParams<false, false> &params, bhcOutputs<false, false> &outputs);
#endif
#if BHC_ENABLE_NX2D
template bool BHC_API
run<true, false>(bhcParams<true, false> &params, bhcOutputs<true, false> &outputs);
#endif
#if BHC_ENABLE_3D
template bool BHC_API
run<true, true>(bhcParams<true, true> &params, bhcOutputs<true, true> &outputs);
#endif

template<bool O3D, bool R3D> bool writeout(
    const bhcParams<O3D, R3D> &params, const bhcOutputs<O3D, R3D> &outputs)
{
    try {
        Stopwatch sw(GetInternal(params));
        sw.tick();
        auto *mo = GetMode<O3D, R3D>(params);
        mo->Writeout(params, outputs);
        sw.tock("writeout");
        delete mo;
    } catch(const std::exception &e) {
        EXTWARN("Exception caught in bhc::writeout(): %s\n", e.what());
        return false;
    }
    return true;
}

#if BHC_ENABLE_2D
template bool BHC_API writeout<false, false>(
    const bhcParams<false, false> &params, const bhcOutputs<false, false> &outputs);
#endif
#if BHC_ENABLE_NX2D
template bool BHC_API writeout<true, false>(
    const bhcParams<true, false> &params, const bhcOutputs<true, false> &outputs);
#endif
#if BHC_ENABLE_3D
template bool BHC_API writeout<true, true>(
    const bhcParams<true, true> &params, const bhcOutputs<true, true> &outputs);
#endif

////////////////////////////////////////////////////////////////////////////////

template<bool O3D, bool R3D> void finalize(
    bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs)
{
    module::ModulesList<O3D, R3D> modules;
    mode::ModesList<O3D, R3D> modes;
    for(auto *m : modules.list()) m->Finalize(params);
    for(auto *m : modes.list()) m->Finalize(params, outputs);

    trackdeallocate(params, params.Bdry);
    trackdeallocate(params, params.bdinfo);
    trackdeallocate(params, params.refl);
    trackdeallocate(params, params.ssp);
    trackdeallocate(params, params.atten);
    trackdeallocate(params, params.Pos);
    trackdeallocate(params, params.Angles);
    trackdeallocate(params, params.freqinfo);
    trackdeallocate(params, params.Beam);
    trackdeallocate(params, params.sbp);
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
