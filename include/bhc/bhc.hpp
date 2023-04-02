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
 *    below.
 * 2. Many data structures have bool flags indicating whether the data is in one
 *    format or another, or whether it has been preprocessed or is dirty. Make
 *    sure you understand any flags present in the structure you are modifying.
 *    For example, params.Pos->SxSyInKm is true if params.Pos->Sx[:] and Sy[:]
 *    are in kilometers and false if they are in meters. If you are writing to
 *    these arrays, you must set the flag to the correct value for your data.
 *    And if you are reading from these arrays, you must read and respect the
 *    current state of the flag. In this case, the data will be converted to
 *    meters (and the flag cleared) when run() or other API calls are made. So,
 *    you may have written the data originally in kilometers, but it may be in
 *    meters now.
 */

/**
 * Reallocate altimetry data to the given size (for 2D this is range, for
 * 3D/Nx2D this is X/Y). After calling this, fill in the coordinates of the
 * boundary points in params.bdinfo->top.bd[:].x. The R/X/Y values need not be
 * uniformly spaced but they must be monotonically increasing.
 * 2D: If params.bdinfo->top.type[1] == 'L', also must fill in params.bdinfo->
 * top.bd[:].hs.{alphaR, betaR, rho, alphaI, betaI}.
 * 3D: Access the boundary array as params.bdinfo->top.bd[ix * NPts.y + iy].
 * The X and Y coordinates must be on a Cartesian grid and filled in into x.x
 * and x.y of each point (even though this duplicates the values many times).
 */
template<bool O3D> void extsetup_altimetry(
    bhcParams<O3D> &params, const IORI2<O3D> &size);
extern template void extsetup_altimetry<false>(
    bhcParams<false> &params, const IORI2<false> &size);
extern template void extsetup_altimetry<true>(
    bhcParams<true> &params, const IORI2<true> &size);
/// See extsetup_altimetry.
template<bool O3D> void extsetup_bathymetry(
    bhcParams<O3D> &params, const IORI2<O3D> &size);
extern template void extsetup_bathymetry<false>(
    bhcParams<false> &params, const IORI2<false> &size);
extern template void extsetup_bathymetry<true>(
    bhcParams<true> &params, const IORI2<true> &size);

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
/// Nx2D or 2D-3D version, see template.
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
 * returns: false if an error occurred, true if no errors.
 */
template<bool O3D, bool R3D> bool writeout(
    const bhcParams<O3D> &params, const bhcOutputs<O3D, R3D> &outputs);

/// 2D version, see template.
extern template BHC_API bool writeout<false, false>(
    const bhcParams<false> &params, const bhcOutputs<false, false> &outputs);
/// Nx2D or 2D-3D version, see template.
extern template BHC_API bool writeout<true, false>(
    const bhcParams<true> &params, const bhcOutputs<true, false> &outputs);
/// 3D version, see template.
extern template BHC_API bool writeout<true, true>(
    const bhcParams<true> &params, const bhcOutputs<true, true> &outputs);

/**
 * Frees memory. You may call run() many times (with changed parameters), you do
 * not have to call setup - run - finalize every time.
 */
template<bool O3D, bool R3D> void finalize(
    bhcParams<O3D> &params, bhcOutputs<O3D, R3D> &outputs);

/// 2D version, see template.
extern template BHC_API void finalize<false, false>(
    bhcParams<false> &params, bhcOutputs<false, false> &outputs);
/// Nx2D or 2D-3D version, see template.
extern template BHC_API void finalize<true, false>(
    bhcParams<true> &params, bhcOutputs<true, false> &outputs);
/// 3D version, see template.
extern template BHC_API void finalize<true, true>(
    bhcParams<true> &params, bhcOutputs<true, true> &outputs);

} // namespace bhc

#ifdef BHC_UNDEF_STD_AFTER
#undef STD
#endif
