#pragma once

#ifndef _BHC_INCLUDED_
#error "This file must be included via #include <bhc/bhc.hpp>!"
#endif

////////////////////////////////////////////////////////////////////////////////
//Select which standard library
////////////////////////////////////////////////////////////////////////////////

#ifdef BHC_BUILD_CUDA
#include "cuda_runtime.h"
#define HOST_DEVICE __host__ __device__
//Requires libcu++ 1.4 or higher, which is included in the normal CUDA Toolkit
//install since somewhere between CUDA 11.3 and 11.5. (I.e. 11.5 and later has
//it, 11.2 and earlier does not, not sure about 11.3 and 11.4. You can check by
//seeing if <cuda_install_dir>/include/cuda/std/complex exists.)
#include <cuda/std/complex>
#include <cuda/std/cfloat>
#define STD cuda::std
#define PROGRAMNAME "bellhopcuda"
#else
#define HOST_DEVICE
#define STD std
#define PROGRAMNAME "bellhopcxx"
#endif

////////////////////////////////////////////////////////////////////////////////
//Shared library setup
////////////////////////////////////////////////////////////////////////////////

#ifdef _WIN32
#ifdef BHC_CMDLINE
#define BHC_API
#elif defined(BHC_EXPORTS)
#define BHC_API __declspec(dllexport)
#else //Users of bellhopcxx / bellhopcuda
#define BHC_API __declspec(dllimport)
#endif
#else
#define BHC_API
#endif
