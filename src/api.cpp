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
#include "module/sbp.hpp"

#include "mode/modemodule.hpp"
#include "mode/ray.hpp"
#include "mode/field.hpp"
#include "mode/tl.hpp"
#include "mode/eigen.hpp"
#include "mode/arr.hpp"

namespace bhc {

namespace module {

template<bool O3D> class ModulesList {
public:
    ModulesList()
    {
        modules.push_back(new Atten<O3D>());
        modules.push_back(new Title<O3D>());
        modules.push_back(new Freq0<O3D>());
        modules.push_back(new NMedia<O3D>());
        modules.push_back(new TopOpt<O3D>());
        modules.push_back(new BoundaryCondTop<O3D>());
        modules.push_back(new SSP<O3D>());
        modules.push_back(new BotOpt<O3D>());
        modules.push_back(new BoundaryCondBot<O3D>());
        modules.push_back(new SxSy<O3D>());
        modules.push_back(new SzRz<O3D>());
        modules.push_back(new RcvrRanges<O3D>());
        modules.push_back(new RcvrBearings<O3D>());
        modules.push_back(new FreqVec<O3D>());
        modules.push_back(new RunType<O3D>());
        modules.push_back(new RayAnglesElevation<O3D>());
        modules.push_back(new RayAnglesBearing<O3D>());
        modules.push_back(new BeamInfo<O3D>());
        modules.push_back(new Altimetry<O3D>());
        modules.push_back(new Bathymetry<O3D>());
        modules.push_back(new BRC<O3D>());
        modules.push_back(new TRC<O3D>());
        modules.push_back(new SBP<O3D>());
    }
    ~ModulesList()
    {
        for(auto *module : modules) delete module;
    }

    const std::vector<ParamsModule<O3D> *> &list() const { return modules; }

private:
    std::vector<ParamsModule<O3D> *> modules;
};

#if BHC_ENABLE_2D
template class ParamsModule<false>;
template class ModulesList<false>;
#endif
#if BHC_ENABLE_NX2D || BHC_ENABLE_3D
template class ParamsModule<true>;
template class ModulesList<true>;
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
template<bool O3D> void setupGPU(const bhcParams<O3D> &params)
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
    const bhcInit &init, bhcParams<O3D> &params, bhcOutputs<O3D, R3D> &outputs)
{
    try {
        params.internal = new bhcInternal(init, O3D, R3D);

        Stopwatch sw(GetInternal(params));
        sw.tick();

        if(GetInternal(params)->maxMemory < 8000000u) {
            EXTERR(
                "%d bytes is an unreasonably small amount of memory to "
                "ask " BHC_PROGRAMNAME " to limit itself to",
                GetInternal(params)->maxMemory);
        }
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

        module::ModulesList<O3D> modules;
        mode::ModesList<O3D, R3D> modes;
        for(auto *m : modules.list()) m->Init(params);
        for(auto *m : modes.list()) m->Init(outputs);

        if(GetInternal(params)->noEnvFil) {
            for(auto *m : modules.list()) {
                m->SetupPre(params);
                m->Default(params);
                m->SetupPost(params);
            }
            for(auto *m : modules.list()) { m->Validate(params); }
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
            for(auto *m : modules.list()) {
                m->Validate(params);
                m->Echo(params);
            }
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
    const bhcInit &init, bhcParams<false> &params, bhcOutputs<false, false> &outputs);
#endif
#if BHC_ENABLE_NX2D
template bool BHC_API setup<true, false>(
    const bhcInit &init, bhcParams<true> &params, bhcOutputs<true, false> &outputs);
#endif
#if BHC_ENABLE_3D
template bool BHC_API setup<true, true>(
    const bhcInit &init, bhcParams<true> &params, bhcOutputs<true, true> &outputs);
#endif

template<bool O3D> bool echo(bhcParams<O3D> &params)
{
    try {
        module::ModulesList<O3D> modules;
        for(auto *m : modules.list()) {
            m->Validate(params);
            m->Echo(params);
        }
    } catch(const std::exception &e) {
        EXTWARN("Exception caught in bhc::echo(): %s\n", e.what());
        return false;
    }

    return true;
}

#if BHC_ENABLE_2D
template BHC_API bool echo<false>(bhcParams<false> &params);
#endif
#if BHC_ENABLE_NX2D || BHC_ENABLE_3D
template BHC_API bool echo<true>(bhcParams<true> &params);
#endif

template<bool O3D> bool writeenv(bhcParams<O3D> &params, const char *FileRoot)
{
    try {
        module::ModulesList<O3D> modules;
        for(auto *m : modules.list()) {
            m->Validate(params);
            m->Preprocess(params);
        }

        GetInternal(params)->FileRoot = FileRoot;
        LDOFile ENVFile;
        ENVFile.setStyle(LDOFile::Style::WRITTEN_BY_HAND);
        ENVFile.open(std::string(FileRoot) + ".env");
        if(!ENVFile.good()) {
            PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;
            PRTFile << "ENVFile = " << FileRoot << ".env\n";
            EXTERR(BHC_PROGRAMNAME
                   " - writeenv: Unable to open the new environmental file");
        }
        for(auto *m : modules.list()) m->Write(params, ENVFile);

    } catch(const std::exception &e) {
        EXTWARN("Exception caught in bhc::writeenv(): %s\n", e.what());
        return false;
    }

    return true;
}

#if BHC_ENABLE_2D
template BHC_API bool writeenv<false>(bhcParams<false> &params, const char *FileRoot);
#endif
#if BHC_ENABLE_NX2D || BHC_ENABLE_3D
template BHC_API bool writeenv<true>(bhcParams<true> &params, const char *FileRoot);
#endif

////////////////////////////////////////////////////////////////////////////////

template<bool O3D, bool R3D> inline mode::ModeModule<O3D, R3D> *GetMode(
    const bhcParams<O3D> &params)
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
        // return nullptr;
    }
}

template<bool O3D, bool R3D> bool run(
    bhcParams<O3D> &params, bhcOutputs<O3D, R3D> &outputs)
{
    try {
        Stopwatch sw(GetInternal(params));

        sw.tick();
        module::ModulesList<O3D> modules;
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
run<false, false>(bhcParams<false> &params, bhcOutputs<false, false> &outputs);
#endif
#if BHC_ENABLE_NX2D
template bool BHC_API
run<true, false>(bhcParams<true> &params, bhcOutputs<true, false> &outputs);
#endif
#if BHC_ENABLE_3D
template bool BHC_API
run<true, true>(bhcParams<true> &params, bhcOutputs<true, true> &outputs);
#endif

template<bool O3D, bool R3D> bool writeout(
    const bhcParams<O3D> &params, const bhcOutputs<O3D, R3D> &outputs,
    const char *FileRoot)
{
    try {
        Stopwatch sw(GetInternal(params));
        sw.tick();
        if(FileRoot != nullptr) { GetInternal(params)->FileRoot = FileRoot; }
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
    const bhcParams<false> &params, const bhcOutputs<false, false> &outputs,
    const char *FileRoot);
#endif
#if BHC_ENABLE_NX2D
template bool BHC_API writeout<true, false>(
    const bhcParams<true> &params, const bhcOutputs<true, false> &outputs,
    const char *FileRoot);
#endif
#if BHC_ENABLE_3D
template bool BHC_API writeout<true, true>(
    const bhcParams<true> &params, const bhcOutputs<true, true> &outputs,
    const char *FileRoot);
#endif

template<bool O3D, bool R3D> bool readout(
    bhcParams<O3D> &params, bhcOutputs<O3D, R3D> &outputs, const char *FileRoot)
{
    try {
        if(FileRoot == nullptr) { FileRoot = GetInternal(params)->FileRoot.c_str(); }
        auto *mo = GetMode<O3D, R3D>(params);
        mo->Readout(params, outputs, FileRoot);
        delete mo;
        module::ModulesList<O3D> modules;
        for(auto *m : modules.list()) m->Validate(params);
    } catch(const std::exception &e) {
        EXTWARN("Exception caught in bhc::readout(): %s\n", e.what());
        return false;
    }
    return true;
}

#if BHC_ENABLE_2D
template BHC_API bool readout<false, false>(
    bhcParams<false> &params, bhcOutputs<false, false> &outputs, const char *FileRoot);
#endif
#if BHC_ENABLE_NX2D
template BHC_API bool readout<true, false>(
    bhcParams<true> &params, bhcOutputs<true, false> &outputs, const char *FileRoot);
#endif
#if BHC_ENABLE_3D
template BHC_API bool readout<true, true>(
    bhcParams<true> &params, bhcOutputs<true, true> &outputs, const char *FileRoot);
#endif

////////////////////////////////////////////////////////////////////////////////

template<bool O3D, bool R3D> void finalize(
    bhcParams<O3D> &params, bhcOutputs<O3D, R3D> &outputs)
{
    module::ModulesList<O3D> modules;
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
template void BHC_API
finalize<false, false>(bhcParams<false> &params, bhcOutputs<false, false> &outputs);
#endif
#if BHC_ENABLE_NX2D
template void BHC_API
finalize<true, false>(bhcParams<true> &params, bhcOutputs<true, false> &outputs);
#endif
#if BHC_ENABLE_3D
template void BHC_API
finalize<true, true>(bhcParams<true> &params, bhcOutputs<true, true> &outputs);
#endif

template<bool O3D> void extsetup_sxsy(bhcParams<O3D> &params, int32_t NSx, int32_t NSy)
{
    module::SxSy<O3D> pm;
    pm.ExtSetup(params, NSx, NSy);
}
#if BHC_ENABLE_2D
template BHC_API void extsetup_sxsy<false>(
    bhcParams<false> &params, int32_t NSx, int32_t NSy);
#endif
#if BHC_ENABLE_3D || BHC_ENABLE_NX2D
template BHC_API void extsetup_sxsy<true>(
    bhcParams<true> &params, int32_t NSx, int32_t NSy);
#endif

template<bool O3D> void extsetup_sz(bhcParams<O3D> &params, int32_t NSz)
{
    module::SzRz<O3D> pm;
    pm.ExtSetupSz(params, NSz);
}
#if BHC_ENABLE_2D
template BHC_API void extsetup_sz<false>(bhcParams<false> &params, int32_t NSz);
#endif
#if BHC_ENABLE_3D || BHC_ENABLE_NX2D
template BHC_API void extsetup_sz<true>(bhcParams<true> &params, int32_t NSz);
#endif

template<bool O3D> void extsetup_rcvrranges(bhcParams<O3D> &params, int32_t NRr)
{
    module::RcvrRanges<O3D> pm;
    pm.ExtSetup(params, NRr);
}
#if BHC_ENABLE_2D
template BHC_API void extsetup_rcvrranges<false>(bhcParams<false> &params, int32_t NRr);
#endif
#if BHC_ENABLE_3D || BHC_ENABLE_NX2D
template BHC_API void extsetup_rcvrranges<true>(bhcParams<true> &params, int32_t NRr);
#endif

template<bool O3D> void extsetup_rcvrdepths(bhcParams<O3D> &params, int32_t NRz)
{
    module::SzRz<O3D> pm;
    pm.ExtSetupRz(params, NRz);
}
#if BHC_ENABLE_2D
template BHC_API void extsetup_rcvrdepths<false>(bhcParams<false> &params, int32_t NRz);
#endif
#if BHC_ENABLE_3D || BHC_ENABLE_NX2D
template BHC_API void extsetup_rcvrdepths<true>(bhcParams<true> &params, int32_t NRz);
#endif

template<bool O3D> void extsetup_rcvrbearings(bhcParams<O3D> &params, int32_t Ntheta)
{
    module::RcvrBearings<O3D> pm;
    pm.ExtSetup(params, Ntheta);
}
#if BHC_ENABLE_2D
template BHC_API void extsetup_rcvrbearings<false>(
    bhcParams<false> &params, int32_t Ntheta);
#endif
#if BHC_ENABLE_3D || BHC_ENABLE_NX2D
template BHC_API void extsetup_rcvrbearings<true>(
    bhcParams<true> &params, int32_t Ntheta);
#endif

template<bool O3D> void extsetup_rayelevations(bhcParams<O3D> &params, int32_t n)
{
    module::RayAnglesElevation<O3D> pm;
    pm.ExtSetup(params, n);
}
#if BHC_ENABLE_2D
template BHC_API void extsetup_rayelevations<false>(bhcParams<false> &params, int32_t n);
#endif
#if BHC_ENABLE_3D || BHC_ENABLE_NX2D
template BHC_API void extsetup_rayelevations<true>(bhcParams<true> &params, int32_t n);
#endif

template<bool O3D> void extsetup_raybearings(bhcParams<O3D> &params, int32_t n)
{
    module::RayAnglesBearing<O3D> pm;
    pm.ExtSetup(params, n);
}
#if BHC_ENABLE_2D
template BHC_API void extsetup_raybearings<false>(bhcParams<false> &params, int32_t n);
#endif
#if BHC_ENABLE_3D || BHC_ENABLE_NX2D
template BHC_API void extsetup_raybearings<true>(bhcParams<true> &params, int32_t n);
#endif

template<bool O3D> void extsetup_sbp(bhcParams<O3D> &params, int32_t NSBPPts)
{
    module::SBP<O3D> pm;
    pm.ExtSetup(params, NSBPPts);
}
#if BHC_ENABLE_2D
template BHC_API void extsetup_sbp<false>(bhcParams<false> &params, int32_t NSBPPts);
#endif
#if BHC_ENABLE_3D || BHC_ENABLE_NX2D
template BHC_API void extsetup_sbp<true>(bhcParams<true> &params, int32_t NSBPPts);
#endif

template<bool O3D> void extsetup_freqvec(bhcParams<O3D> &params, int32_t Nfreq)
{
    module::FreqVec<O3D> pm;
    pm.ExtSetup(params, Nfreq);
}
#if BHC_ENABLE_2D
template BHC_API void extsetup_freqvec<false>(bhcParams<false> &params, int32_t Nfreq);
#endif
#if BHC_ENABLE_3D || BHC_ENABLE_NX2D
template BHC_API void extsetup_freqvec<true>(bhcParams<true> &params, int32_t Nfreq);
#endif

template<bool O3D> void extsetup_altimetry(bhcParams<O3D> &params, const IORI2<O3D> &NPts)
{
    module::Altimetry<O3D> pm;
    pm.ExtSetup(params, NPts);
}
#if BHC_ENABLE_2D
template BHC_API void extsetup_altimetry<false>(
    bhcParams<false> &params, const IORI2<false> &NPts);
#endif
#if BHC_ENABLE_3D || BHC_ENABLE_NX2D
template BHC_API void extsetup_altimetry<true>(
    bhcParams<true> &params, const IORI2<true> &NPts);
#endif

template<bool O3D> void extsetup_bathymetry(
    bhcParams<O3D> &params, const IORI2<O3D> &NPts)
{
    module::Bathymetry<O3D> pm;
    pm.ExtSetup(params, NPts);
}
#if BHC_ENABLE_2D
template BHC_API void extsetup_bathymetry<false>(
    bhcParams<false> &params, const IORI2<false> &NPts);
#endif
#if BHC_ENABLE_3D || BHC_ENABLE_NX2D
template BHC_API void extsetup_bathymetry<true>(
    bhcParams<true> &params, const IORI2<true> &NPts);
#endif

template<bool O3D> void extsetup_trc(bhcParams<O3D> &params, int32_t NPts)
{
    module::TRC<O3D> pm;
    pm.ExtSetup(params, NPts);
}
#if BHC_ENABLE_2D
template BHC_API void extsetup_trc<false>(bhcParams<false> &params, int32_t NPts);
#endif
#if BHC_ENABLE_3D || BHC_ENABLE_NX2D
template BHC_API void extsetup_trc<true>(bhcParams<true> &params, int32_t NPts);
#endif

template<bool O3D> void extsetup_brc(bhcParams<O3D> &params, int32_t NPts)
{
    module::BRC<O3D> pm;
    pm.ExtSetup(params, NPts);
}
#if BHC_ENABLE_2D
template BHC_API void extsetup_brc<false>(bhcParams<false> &params, int32_t NPts);
#endif
#if BHC_ENABLE_3D || BHC_ENABLE_NX2D
template BHC_API void extsetup_brc<true>(bhcParams<true> &params, int32_t NPts);
#endif

#if BHC_ENABLE_2D
extern BHC_API void extsetup_ssp_quad(bhcParams<false> &params, int32_t NPts, int32_t Nr)
{
    module::SSP<false> pm;
    pm.ExtSetup(params, NPts, Nr, 0);
}
#endif
#if BHC_ENABLE_3D || BHC_ENABLE_NX2D
extern BHC_API void extsetup_ssp_hexahedral(
    bhcParams<true> &params, int32_t Nx, int32_t Ny, int32_t Nz)
{
    module::SSP<true> pm;
    pm.ExtSetup(params, Nx, Ny, Nz);
}
#endif

} // namespace bhc
