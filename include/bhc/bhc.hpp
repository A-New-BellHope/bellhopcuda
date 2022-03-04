#pragma once

#define _USE_MATH_DEFINES 1 //must be before anything which includes math.h
#include <math.h>

#define _BHC_INCLUDED_ 1
#include "platform.hpp"
#include "real.hpp"
#include "cpx.hpp"
#include "structs.hpp"
#undef _BHC_INCLUDED_

#include <iostream>
#include <fstream>

/**
 * Optional PRTFile initialization. You can also use an ostringstream for
 * PRTFile but then you have to set it up yourself.
 */
static inline void OpenPRTFile(std::string FileRoot, std::ofstream &PRTFile){
    PRTFile.open(FileRoot + ".prt");
    if(!PRTFile.good()){
        std::cout << "Could not open print file: " << FileRoot << ".prt\n";
        std::abort();
    }
    PRTFile << std::unitbuf;
}

/**
 * Main BELLHOP setup from an environment file. Call this to create and
 * initialize the params. You may modify the params after calling this and
 * before calling run().
 * 
 * FileRoot: Relative path to environment file, without the .env extension. E.g. 
 * path/to/MunkB_ray_rot (where path/to/MunkB_ray_rot.env and also path/to/
 * MunkB_ray_rot.ssp, path/to/MunkB_ray_rot.bty, etc. exist).
 * 
 * PRTFile: Where BELLHOP debug / informational output is sent. This may be an
 * ofstream or an ostringstream.
 * 
 * params, outputs: Just create uninitialized structs and pass them in to be
 * initialized.
 */
BHC_API void setup(std::string FileRoot, std::ostream &PRTFile, bhcParams &params,
    bhcOutputs &outputs);
    
/**
 * Runs the selected run type and places the results in the appropriate struct
 * within outputs. outputs need not be initialized prior to the call.
 */
BHC_API void run(std::ostream &PRTFile, const bhcParams &params, bhcOutputs &outputs,
    bool singlethread);

/**
 * Frees memory. You may call run() many times, you do not have to call setup
 * - run - finalize every time.
 */
//BHC_API void finalize(bhcParams &params, bhcOutputs &outputs); TODO
