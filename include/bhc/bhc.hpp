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
 */
BHC_API bool setup(const char *FileRoot, void (*outputCallback)(const char *message),
    bhcParams &params, bhcOutputs &outputs);
    
/**
 * Runs the selected run type and places the results in the appropriate struct
 * within outputs. outputs need not be initialized prior to the call.
 * 
 * An env file should usually be read directly first, by calling setup then 
 * run-finalize. But, before calling finalize, you may edit parameters and rerun 
 * with an expected pattern
 *     setup-run-change params-run-change params...-finalize.
 * TODO: Only a few parameters can be updated, notably sources and 1D SSP - JS, 8/25/2022.
 * 
 * returns: false on fatal errors, true otherwise. If a fatal error occurs,
 * must call finalize() and setup() again before continuing to use the library.
 */
BHC_API bool run(bhcParams &params, bhcOutputs &outputs, bool singlethread);

/**
 * Frees memory. You may call run() many times, you do not have to call setup
 * - run - finalize every time.
 */
BHC_API void finalize(bhcParams &params, bhcOutputs &outputs);

}

#ifdef BHC_UNDEF_STD_AFTER
#undef STD
#endif
