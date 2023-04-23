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
#pragma once

#define _USE_MATH_DEFINES 1 // must be before anything which includes math.h
#include <math.h>

#define _BHC_INCLUDED_ 1
#include "platform.hpp"
#include "math.hpp"
#include "structs.hpp"
#undef _BHC_INCLUDED_

namespace bhc {

/**
 * NOTE: If you are on Windows and writing a program which will link to the
 * bellhopcxx / bellhopcuda DLL, you must define BHC_DLL_IMPORT before including
 * this header.
 *
 * Main BELLHOP setup from an environment file. Call this to create and
 * initialize the params. You may modify the params after calling this and
 * before calling run().
 *
 * You may use "multiple instances" of bellhopcxx / bellhopcuda within the same
 * process by calling this (and the other functions below) with different params
 * and outputs; there are no global variables in the library.
 *
 * init: Initialization parameters. See the documentation of each of the members
 * of the struct for more info.
 *
 * params, outputs: Just create uninitialized structs and pass them in to be
 * initialized. You may modify params after setup.
 *
 * returns: false if an error occurred, true if no errors.
 *
 * O3D stands for "ocean 3D" and R3D stands for "ray(s) 3D".
 * O3D=false, R3D=false: 2D mode
 * O3D=true, R3D=false: Nx2D mode
 * O3D=true, R3D=true: 3D mode
 */
template<bool O3D, bool R3D> bool setup(
    const bhcInit &init, bhcParams<O3D> &params, bhcOutputs<O3D, R3D> &outputs);

/// 2D version, see template.
extern template BHC_API bool setup<false, false>(
    const bhcInit &init, bhcParams<false> &params, bhcOutputs<false, false> &outputs);
/// Nx2D version, see template.
extern template BHC_API bool setup<true, false>(
    const bhcInit &init, bhcParams<true> &params, bhcOutputs<true, false> &outputs);
/// 3D version, see template.
extern template BHC_API bool setup<true, true>(
    const bhcInit &init, bhcParams<true> &params, bhcOutputs<true, true> &outputs);

/*
 * You can generally modify params as desired before run. There are two main
 * restrictions on this.
 * 1. You must not allocate or deallocate any data structures / arrays within
 *    params. If you need to change their size, you must use the functions
 *    below to reallocate the structures to the desired size, and then fill in
 *    the data into the arrays.
 * 2. Many data structures have bool flags indicating whether the data is in one
 *    format or another, or whether it has been preprocessed or is dirty. Make
 *    sure you understand any flags present in the structure you are modifying.
 *    For example, params.Pos->SxSyInKm is true if params.Pos->Sx[:] and Sy[:]
 *    are in kilometers and false if they are in meters. If you are writing to
 *    these arrays, you must set the flag to the correct value for your data
 *    (don't assume its current state). And if you are reading from these
 *    arrays, you must read and respect the current state of the flag. In this
 *    case, the data will be converted to meters (and the flag cleared) when
 *    run() or other API calls are made. So, you may have written the data
 *    originally in kilometers, but it may be in meters now.
 * Also, most arrays must be monotonically increasing along relevant axes (e.g.
 * bathymetry X and Y values must be monotonic but Z values can be arbitrary).
 */

/**
 * Reallocate the source X and Y positions to the given size. In 2D mode, there
 * must be exactly one source in each direction at coordinate 0.0; it is set up
 * this way by default so you never need to call this function. Set
 * params.Pos->SxSyInKm depending on whether you write the data in kilometers or
 * meters.
 */
template<bool O3D> void extsetup_sxsy(bhcParams<O3D> &params, int32_t NSx, int32_t NSy);
extern template BHC_API void extsetup_sxsy<false>(
    bhcParams<false> &params, int32_t NSx, int32_t NSy);
extern template BHC_API void extsetup_sxsy<true>(
    bhcParams<true> &params, int32_t NSx, int32_t NSy);
/**
 * Reallocate the source Z positions (depths) to the given size. Depth values
 * are always specified in meters.
 */
template<bool O3D> void extsetup_sz(bhcParams<O3D> &params, int32_t NSz);
extern template BHC_API void extsetup_sz<false>(bhcParams<false> &params, int32_t NSz);
extern template BHC_API void extsetup_sz<true>(bhcParams<true> &params, int32_t NSz);
/**
 * Reallocate the receiver ranges to the given size. Set params.Pos->RrInKm
 * depending on whether you enter the data in kilometers or meters.
 */
template<bool O3D> void extsetup_rcvrranges(bhcParams<O3D> &params, int32_t NRr);
extern template BHC_API void extsetup_rcvrranges<false>(
    bhcParams<false> &params, int32_t NRr);
extern template BHC_API void extsetup_rcvrranges<true>(
    bhcParams<true> &params, int32_t NRr);
/**
 * Reallocate the receiver Z positions (depths) to the given size. Depth values
 * are always specified in meters.
 */
template<bool O3D> void extsetup_rcvrdepths(bhcParams<O3D> &params, int32_t NRz);
extern template BHC_API void extsetup_rcvrdepths<false>(
    bhcParams<false> &params, int32_t NRz);
extern template BHC_API void extsetup_rcvrdepths<true>(
    bhcParams<true> &params, int32_t NRz);
/**
 * Reallocate the receiver bearing angles to the given size. In 2D mode, there
 * must be exactly one bearing angle of 0.0; it is set up this way by default so
 * you never need to call this function. Receiver bearings are always in degrees.
 */
template<bool O3D> void extsetup_rcvrbearings(bhcParams<O3D> &params, int32_t Ntheta);
extern template BHC_API void extsetup_rcvrbearings<false>(
    bhcParams<false> &params, int32_t Ntheta);
extern template BHC_API void extsetup_rcvrbearings<true>(
    bhcParams<true> &params, int32_t Ntheta);
/**
 * Reallocate the ray elevation angles to the given size. Set
 * params.Angles->alpha.inDegrees depending on whether you enter the data in
 * degrees or radians.
 */
template<bool O3D> void extsetup_rayelevations(bhcParams<O3D> &params, int32_t n);
extern template BHC_API void extsetup_rayelevations<false>(
    bhcParams<false> &params, int32_t n);
extern template BHC_API void extsetup_rayelevations<true>(
    bhcParams<true> &params, int32_t n);
/**
 * Reallocate the ray bearing angles to the given size. In 2D mode, there must
 * be exactly one bearing angle of 0.0; it is set up this way by default so you
 * never need to call this function. Note that in Nx2D mode, these angles are
 * replaced with the receiver bearing angles. Set params.Angles->beta.inDegrees
 * depending on whether you enter the data in degrees or radians.
 */
template<bool O3D> void extsetup_raybearings(bhcParams<O3D> &params, int32_t n);
extern template BHC_API void extsetup_raybearings<false>(
    bhcParams<false> &params, int32_t n);
extern template BHC_API void extsetup_raybearings<true>(
    bhcParams<true> &params, int32_t n);
/**
 * Reallocate the source beam pattern to the given size. After calling this,
 * write angles (monotonic, -180.0 to 180.0) and levels into params.sbp->SrcBmPat
 * with element 0 being the angle, element 1 being the level, and so on.
 * Set params.sbp->SBPIndB based on whether the levels are in dB or linear
 * amplitude.
 */
template<bool O3D> void extsetup_sbp(bhcParams<O3D> &params, int32_t NSBPPts);
extern template BHC_API void extsetup_sbp<false>(
    bhcParams<false> &params, int32_t NSBPPts);
extern template BHC_API void extsetup_sbp<true>(bhcParams<true> &params, int32_t NSBPPts);
/**
 * Reallocate the broadband frequency vector to the given size. FreqVec
 * (broadband mode) is not properly supported in BELLHOP(3D)--these values have
 * no impact on the physics model--but this feature cannot be removed from
 * bellhopcxx / bellhopcuda because the frequency vector is written to the
 * SHDFile.
 */
template<bool O3D> void extsetup_freqvec(bhcParams<O3D> &params, int32_t Nfreq);
extern template BHC_API void extsetup_freqvec<false>(
    bhcParams<false> &params, int32_t Nfreq);
extern template BHC_API void extsetup_freqvec<true>(
    bhcParams<true> &params, int32_t Nfreq);
/**
 * Reallocate altimetry data to the given size (for 2D this is range, for
 * 3D/Nx2D this is X/Y). After calling this, fill in the coordinates of the
 * boundary points in params.bdinfo->top.bd[:].x. The R/X/Y values need not be
 * uniformly spaced but they must be monotonically increasing.
 *
 * 2D: If params.bdinfo->top.type[1] == 'L', also must fill in params.bdinfo->
 * top.bd[:].hs.{alphaR, betaR, rho, alphaI, betaI}.
 *
 * 3D: Access the boundary array as params.bdinfo->top.bd[ix * NPts.y + iy].
 * The X and Y coordinates must be on a Cartesian grid and filled in into x.x
 * and x.y of each point (even though this duplicates the values many times).
 *
 * This function sets params.bdinfo->top.dirty, but if you change the altimetry
 * data later (not immediately after calling this function), you'll need to set
 * the dirty flag each time it is changed.
 */
template<bool O3D> void extsetup_altimetry(
    bhcParams<O3D> &params, const IORI2<O3D> &NPts);
extern template BHC_API void extsetup_altimetry<false>(
    bhcParams<false> &params, const IORI2<false> &NPts);
extern template BHC_API void extsetup_altimetry<true>(
    bhcParams<true> &params, const IORI2<true> &NPts);
/// See extsetup_altimetry.
template<bool O3D> void extsetup_bathymetry(
    bhcParams<O3D> &params, const IORI2<O3D> &NPts);
extern template BHC_API void extsetup_bathymetry<false>(
    bhcParams<false> &params, const IORI2<false> &NPts);
extern template BHC_API void extsetup_bathymetry<true>(
    bhcParams<true> &params, const IORI2<true> &NPts);
/**
 * Reallocate the top reflection coefficients to the given size. This also sets
 * params.Bdry->Top.hs.Opt[1] to 'F' (file) to use the reflection coefficients.
 * If you don't want to use them anymore, set that flag to something else like
 * 'V' (vacuum). Also set params.refl->top.inDegrees to the appropriate value.
 */
template<bool O3D> void extsetup_trc(bhcParams<O3D> &params, int32_t NPts);
extern template BHC_API void extsetup_trc<false>(bhcParams<false> &params, int32_t NPts);
extern template BHC_API void extsetup_trc<true>(bhcParams<true> &params, int32_t NPts);
/**
 * Reallocate the bottom reflection coefficients to the given size. This also sets
 * params.Bdry->Bot.hs.Opt[0] to 'F' (file) to use the reflection coefficients.
 * If you don't want to use them anymore, set that flag to something else like
 * 'R' (rigid). Also set params.refl->bot.inDegrees to the appropriate value.
 */
template<bool O3D> void extsetup_brc(bhcParams<O3D> &params, int32_t NPts);
extern template BHC_API void extsetup_brc<false>(bhcParams<false> &params, int32_t NPts);
extern template BHC_API void extsetup_brc<true>(bhcParams<true> &params, int32_t NPts);
/**
 * Set up and/or reallocate the SSP for quad mode (2D only). NPts is the number
 * of depths. Fill in params.ssp->z, params.ssp->Seg.r, and params.ssp->cMat[z *
 * Nr + r]. Set params.ssp->rangeInKm depending on whether you are specifying
 * ranges in kilometers or meters.
 *
 * After writing your SSP, make sure params.Bdry->Top.hs.Depth (nominal surface
 * depth, normally zero) is equal to ssp->z[0], and params.Bdry->Bot.hs.Depth
 * (nominal ocean bottom depth) is equal to ssp->z[NPts-1].
 *
 * This function sets params.ssp->dirty, but if you change the SSP data later
 * (not immediately after calling this function), you'll need to set the dirty
 * flag each time it is changed.
 *
 * To set the SSP to a mode other than quad or hexahedral, an extsetup call is
 * not needed. Just
 * - set params.ssp->Type to the correct letter
 * - set params.ssp->NPts and fill in params->ssp.z[:], .alphaR, .betaR, .rho,
 *   .alphaI, and .betaI
 * - make sure the surface and bottom depths are set correctly as described
 *   above
 * - set the dirty flag
 */
extern BHC_API void extsetup_ssp_quad(bhcParams<false> &params, int32_t NPts, int32_t Nr);
/**
 * Set up and/or reallocate the SSP for hexahedral mode (3D/Nx2D only). Fill in
 * the X, Y, and Z coordinates (all monotonically increasing, but not
 * necessarily uniformly sampled) in params.ssp->Seg.x, .y, .z, and the speeds
 * in params.ssp->cMat[(x*Ny+y)*Nz+z]. Set params.ssp->rangeInKm depending on
 * whether you are specifying ranges in kilometers or meters.
 *
 * After writing your SSP, make sure params.Bdry->Top.hs.Depth (nominal surface
 * depth, normally zero) is equal to ssp->Seg.z[0], and
 * params.Bdry->Bot.hs.Depth (nominal ocean bottom depth) is equal to
 * ssp->Seg.z[Nz-1].
 *
 * This function sets params.ssp->dirty, but if you change the SSP data later
 * (not immediately after calling this function), you'll need to set the dirty
 * flag each time it is changed.
 *
 * To set the SSP to a mode other than quad or hexahedral, an extsetup call is
 * not needed; see the doc for extsetup_ssp_quad() for more info.
 */
extern BHC_API void extsetup_ssp_hexahedral(
    bhcParams<true> &params, int32_t Nx, int32_t Ny, int32_t Nz);

/**
 * Validates the state of params and writes a summary of the state to the
 * PRTFile or callback. This is done automatically as part of setup() if you
 * started from an environment file, but if you started from defaults and then
 * wrote your own data in, you might want to check whether that data is correct.
 * (The validation step is also performed as part of run().)
 */
template<bool O3D> bool echo(bhcParams<O3D> &params);

/// 2D version, see template.
extern template BHC_API bool echo<false>(bhcParams<false> &params);
/// 3D or Nx2D version, see template.
extern template BHC_API bool echo<true>(bhcParams<true> &params);

/**
 * Runs the selected run type and places the results in the appropriate struct
 * within outputs.
 *
 * returns: false if an error occurred, true if no errors.
 */
template<bool O3D, bool R3D> bool run(
    bhcParams<O3D> &params, bhcOutputs<O3D, R3D> &outputs);

/// 2D version, see template.
extern template BHC_API bool run<false, false>(
    bhcParams<false> &params, bhcOutputs<false, false> &outputs);
/// Nx2D version, see template.
extern template BHC_API bool run<true, false>(
    bhcParams<true> &params, bhcOutputs<true, false> &outputs);
/// 3D version, see template.
extern template BHC_API bool run<true, true>(
    bhcParams<true> &params, bhcOutputs<true, true> &outputs);

/**
 * Write results for the past run to BELLHOP-formatted files, i.e. a ray file,
 * a shade file, or an arrivals file. If you only want to use the results in
 * memory, there is no need to call this.
 *
 * You can pass nullptr for FileRoot to write to output files relative to the
 * same FileRoot that the environment file was originally loaded from.
 *
 * returns: false if an error occurred, true if no errors.
 */
template<bool O3D, bool R3D> bool writeout(
    const bhcParams<O3D> &params, const bhcOutputs<O3D, R3D> &outputs,
    const char *FileRoot);

/// 2D version, see template.
extern template BHC_API bool writeout<false, false>(
    const bhcParams<false> &params, const bhcOutputs<false, false> &outputs,
    const char *FileRoot);
/// Nx2D version, see template.
extern template BHC_API bool writeout<true, false>(
    const bhcParams<true> &params, const bhcOutputs<true, false> &outputs,
    const char *FileRoot);
/// 3D version, see template.
extern template BHC_API bool writeout<true, true>(
    const bhcParams<true> &params, const bhcOutputs<true, true> &outputs,
    const char *FileRoot);

/**
 * Read saved results from a past run (a ray file, TL / shade file, or arrivals
 * file) to memory (the outputs struct). params should have already been
 * initialized with the same or a similar environment file (and other input data
 * files) which produced these outputs. If basic parameters in the output file
 * (e.g number of sources) contradict the state in params, params will be
 * updated to match the output file. However, if there are major
 * incompatibilities such as dimensionality or run type, errors will be thrown.
 *
 * You can pass nullptr for FileRoot to read from output files relative to the
 * same FileRoot that the environment file was originally loaded from.
 *
 * returns: false if an error occurred, true if no errors.
 */
template<bool O3D, bool R3D> bool readout(
    bhcParams<O3D> &params, bhcOutputs<O3D, R3D> &outputs, const char *FileRoot);

/// 2D version, see template.
extern template BHC_API bool readout<false, false>(
    bhcParams<false> &params, bhcOutputs<false, false> &outputs, const char *FileRoot);
/// Nx2D version, see template.
extern template BHC_API bool readout<true, false>(
    bhcParams<true> &params, bhcOutputs<true, false> &outputs, const char *FileRoot);
/// 3D version, see template.
extern template BHC_API bool readout<true, true>(
    bhcParams<true> &params, bhcOutputs<true, true> &outputs, const char *FileRoot);

/**
 * Write the current params state to an environment file and any other "input"
 * file types (e.g. SSP, bathymetry, etc.), so that the state can be loaded
 * later.
 *
 * returns: false if an error occurred, true if no errors.
 */
template<bool O3D> bool writeenv(bhcParams<O3D> &params, const char *FileRoot);

/// 2D version, see template.
extern template BHC_API bool writeenv<false>(
    bhcParams<false> &params, const char *FileRoot);
/// 3D or Nx2D version, see template.
extern template BHC_API bool writeenv<true>(
    bhcParams<true> &params, const char *FileRoot);

/**
 * Frees memory. You may call run() many times (with changed parameters), you do
 * not have to call setup - run - finalize every time.
 */
template<bool O3D, bool R3D> void finalize(
    bhcParams<O3D> &params, bhcOutputs<O3D, R3D> &outputs);

/// 2D version, see template.
extern template BHC_API void finalize<false, false>(
    bhcParams<false> &params, bhcOutputs<false, false> &outputs);
/// Nx2D version, see template.
extern template BHC_API void finalize<true, false>(
    bhcParams<true> &params, bhcOutputs<true, false> &outputs);
/// 3D version, see template.
extern template BHC_API void finalize<true, true>(
    bhcParams<true> &params, bhcOutputs<true, true> &outputs);

} // namespace bhc

#ifdef BHC_UNDEF_STD_AFTER
#undef STD
#endif
