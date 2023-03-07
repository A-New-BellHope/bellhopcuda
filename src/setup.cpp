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

constexpr bool Init_Inline = false;

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
        HSInfo RecycledHS; // Values only initialized once--reused from top to ssp, and
                           // ssp to bot

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
        params.Pos->iSz                = nullptr;
        params.Pos->iRz                = nullptr;
        params.Pos->Sx                 = nullptr;
        params.Pos->Sy                 = nullptr;
        params.Pos->Sz                 = nullptr;
        params.Pos->Rr                 = nullptr;
        params.Pos->Rz                 = nullptr;
        params.Pos->ws                 = nullptr;
        params.Pos->wr                 = nullptr;
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
        RecycledHS.alphaR                = FL(1500.0);
        RecycledHS.betaR                 = FL(0.0);
        RecycledHS.alphaI                = FL(0.0);
        RecycledHS.betaI                 = FL(0.0);
        RecycledHS.rho                   = FL(1.0);

        if constexpr(!O3D && !R3D && Init_Inline) {
            // NPts, Sigma not used by BELLHOP
            std::string TempTitle = BHC_PROGRAMNAME
                "- Calibration case with envfil passed as parameters";
            size_t l = bhc::min(sizeof(params.Title) - 1, TempTitle.size());
            memcpy(params.Title, TempTitle.c_str(), l);
            params.Title[l]        = 0;
            params.freqinfo->freq0 = FL(250.0);
            // NMedia variable is not used by BELLHOP

            // *** Boundary information (type of boundary condition and, if a halfspace,
            // then halfspace info)

            memcpy(params.ssp->AttenUnit, "W", 2); // LP: not a typo--one character string
                                                   // assigned to two
            params.Bdry->Top.hs.bc    = 'V';
            params.Bdry->Top.hs.Depth = FL(0.0);
            params.Bdry->Bot.hs.Depth = FL(100.0);
            memcpy(params.Bdry->Bot.hs.Opt, "A_", 2);
            params.Bdry->Bot.hs.bc = 'A';
            // compressional wave speed
            params.Bdry->Bot.hs.cP
                = crci(params, RL(1.0e20), RL(1590.0), RL(0.5), params.ssp->AttenUnit);
            // shear         wave speed
            params.Bdry->Bot.hs.cS
                = crci(params, RL(1.0e20), RL(0.0), RL(0.0), params.ssp->AttenUnit);
            params.Bdry->Bot.hs.rho = FL(1.2);

            // *** sound speed in the water column ***

            params.ssp->Type  = 'C'; // interpolation method for SSP
            params.ssp->NPts  = 2;   // number of SSP points
            params.ssp->z[0]  = FL(0.0);
            params.ssp->z[1]  = FL(100.0);
            params.ssp->c[0]  = FL(1500.0);
            params.ssp->c[1]  = FL(1500.0);
            params.ssp->cz[0] = FL(0.0);
            params.ssp->cz[1] = FL(0.0); // user should really not have to supply this ...

            // *** source and receiver positions ***

            params.Pos->NSz = 1;
            params.Pos->NRz = 100;
            params.Pos->NRr = 500;

            trackallocate(
                params, "unused inline arrays", params.Pos->Sz, params.Pos->NSz);
            trackallocate(
                params, "unused inline arrays", params.Pos->ws, params.Pos->NSz);
            trackallocate(
                params, "unused inline arrays", params.Pos->iSz, params.Pos->NSz);
            trackallocate(
                params, "unused inline arrays", params.Pos->Rz, params.Pos->NRz);
            trackallocate(
                params, "unused inline arrays", params.Pos->wr, params.Pos->NRz);
            trackallocate(
                params, "unused inline arrays", params.Pos->iRz, params.Pos->NRz);
            trackallocate(
                params, "unused inline arrays", params.Pos->Rr, params.Pos->NRr);

            params.Pos->Sz[0] = FL(50.0);
            // params.Pos->Rz = {0, 50, 100};
            // params.Pos->Rr = {1000, 2000, 3000, 4000, 5000}; // meters !!!
            for(int32_t jj = 0; jj < params.Pos->NRz; ++jj) params.Pos->Rz[jj] = jj;
            for(int32_t jj = 0; jj < params.Pos->NRr; ++jj)
                params.Pos->Rr[jj] = FL(10.0) * jj;

            memcpy(params.Beam->RunType, "C      ", 7);
            memcpy(params.Beam->Type, "G   ", 4);
            params.Beam->deltas = FL(0.0);
            params.Beam->Box.y  = FL(101.0);
            params.Beam->Box.x  = FL(5100.0); // meters

            params.Angles->alpha.n = 1789;
            // params.Angles->alpha.angles = {-80, -70, -60, -50, -40, -30, -20, -10, 0,
            // 10, 20, 30, 40, 50, 60, 70, 80}; // -89 89
            for(int32_t jj = 0; jj < params.Angles->alpha.n; ++jj) {
                params.Angles->alpha.angles[jj] = (FL(180.0) / params.Angles->alpha.n)
                        * (real)jj
                    - FL(90.0);
            }

            // *** altimetry ***

            trackallocate(params, "unused inline arrays", params.bdinfo->top.bd, 2);
            params.bdinfo->top.bd[0].x = vec2(-bdry_big<false>::value(), RL(0.0));
            params.bdinfo->top.bd[1].x = vec2(bdry_big<false>::value(), RL(0.0));

            ComputeBdryTangentNormal<O3D>(&params.bdinfo->top, true);

            // *** bathymetry ***

            trackallocate(params, "unused inline arrays", params.bdinfo->bot.bd, 2);
            params.bdinfo->bot.bd[0].x = vec2(-bdry_big<false>::value(), RL(5000.0));
            params.bdinfo->bot.bd[1].x = vec2(bdry_big<false>::value(), RL(5000.0));

            ComputeBdryTangentNormal<O3D>(&params.bdinfo->bot, false);

            trackallocate(params, "unused inline arrays", params.refl->bot.r, 1);
            trackallocate(params, "unused inline arrays", params.refl->top.r, 1);
            params.refl->bot.NPts = params.refl->top.NPts = 1;

            // *** Source Beam Pattern ***
            params.beaminfo->NSBPPts = 2;
            trackallocate(
                params, "unused inline arrays", params.beaminfo->SrcBmPat, 2 * 2);
            params.beaminfo->SrcBmPat[0 * 2 + 0] = FL(-180.0);
            params.beaminfo->SrcBmPat[0 * 2 + 1] = FL(0.0);
            params.beaminfo->SrcBmPat[1 * 2 + 0] = FL(180.0);
            params.beaminfo->SrcBmPat[1 * 2 + 1] = FL(0.0);
            for(int32_t i = 0; i < 2; ++i)
                params.beaminfo->SrcBmPat[i * 2 + 1] = STD::pow(
                    FL(10.0), params.beaminfo->SrcBmPat[i * 2 + 1] / FL(20.0)); // convert
                                                                                // dB to
                                                                                // linear
                                                                                // scale
                                                                                // !!!
        } else {
            ReadEnvironment<O3D, R3D>(params, RecycledHS);
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
