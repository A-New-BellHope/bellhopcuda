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

#define _USE_MATH_DEFINES 1 //must be before anything which includes math.h
#include <math.h>

#define _BHC_INCLUDED_ 1
#include "platform.hpp"
#include "math.hpp"
#include "structs.hpp"
#undef _BHC_INCLUDED_

namespace bhc {

/**
 * Main BELLHOP setup from an environment file. Call this to create and
 * initialize the params. You may modify the params after calling this and
 * before calling run().
 * 
 * FileRoot: Relative path to environment file, without the .env extension. E.g. 
 * path/to/MunkB_ray_rot (where path/to/MunkB_ray_rot.env and also path/to/
 * MunkB_ray_rot.ssp, path/to/MunkB_ray_rot.bty, etc. exist).
 * 
 * outputCallback: Callback called by setup/run code which will be called for
 * messages (e.g. debug output, error messages). If nullptr is passed, will
 * open a PRTFile (<FileRoot>.prt) and put the messages in there.
 * 
 * params, outputs: Just create uninitialized structs and pass them in to be
 * initialized.
 * 
 * returns: false on fatal errors, true otherwise. If a fatal error occurs,
 * must call finalize() and setup() again before continuing to use the library.
 * 
 * O3D stands for "ocean 3D" and R3D stands for "ray(s) 3D".
 * O3D=false, R3D=false: 2D mode
 * O3D=true, R3D=false: Nx2D mode
 * O3D=true, R3D=true: 3D mode
 */
template<bool O3D, bool R3D> bool setup(
    const char *FileRoot, void (*outputCallback)(const char *message),
    bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs);

/// 2D version, see template.
BHC_API extern template bool setup<false, false>(
    const char *FileRoot, void (*outputCallback)(const char *message),
    bhcParams<false, false> &params, bhcOutputs<false, false> &outputs);
/// Nx2D or 2D-3D version, see template.
BHC_API extern template bool setup<true, false>(
    const char *FileRoot, void (*outputCallback)(const char *message),
    bhcParams<true, false> &params, bhcOutputs<true, false> &outputs);
/// 3D version, see template.
BHC_API extern template bool setup<true, true>(
    const char *FileRoot, void (*outputCallback)(const char *message),
    bhcParams<true, true> &params, bhcOutputs<true, true> &outputs);
    
/**
 * Runs the selected run type and places the results in the appropriate struct
 * within outputs. outputs need not be initialized prior to the call.
 * 
 * returns: false on fatal errors, true otherwise. If a fatal error occurs,
 * must call finalize() and setup() again before continuing to use the library.
 * 
 * Don't call this from multiple threads at the same time (e.g. with different
 * parameters); there is only one static copy of the functionality for
 * synchronizing the different threads launched by these functions, so this will
 * not work correctly from multiple calls simultaneously. Don't do this even
 * if singlethread is set.
 */
template<bool O3D, bool R3D> bool run(
    const bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs, bool singlethread);

/// 2D version, see template.
BHC_API extern template bool run<false, false>(
    const bhcParams<false, false> &params, bhcOutputs<false, false> &outputs, bool singlethread);
/// Nx2D or 2D-3D version, see template.
BHC_API extern template bool run<true, false>(
    const bhcParams<true, false> &params, bhcOutputs<true, false> &outputs, bool singlethread);
/// 3D version, see template.
BHC_API extern template bool run<true, true>(
    const bhcParams<true, true> &params, bhcOutputs<true, true> &outputs, bool singlethread); 

/**
 * Frees memory. You may call run() many times, you do not have to call setup
 * - run - finalize every time.
 */
template<bool O3D, bool R3D> void finalize(
    bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs);

/// 2D version, see template.
BHC_API extern template void finalize<false, false>(
    bhcParams<false, false> &params, bhcOutputs<false, false> &outputs);
/// Nx2D or 2D-3D version, see template.
BHC_API extern template void finalize<true, false>(
    bhcParams<true, false> &params, bhcOutputs<true, false> &outputs);
/// 3D version, see template.
BHC_API extern template void finalize<true, true>(
    bhcParams<true, true> &params, bhcOutputs<true, true> &outputs);

}

#ifdef BHC_UNDEF_STD_AFTER
#undef STD
#endif
