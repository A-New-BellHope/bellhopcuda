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

////////////////////////////////////////////////////////////////////////////////
// General headers
////////////////////////////////////////////////////////////////////////////////

#define _USE_MATH_DEFINES 1 // must be before anything which includes math.h
#include <math.h>

#ifdef _MSC_VER
#define NOMINMAX
#include <windows.h>
#endif

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstring>
#include <cinttypes>
#include <cstdarg>
#include <chrono>
#include <thread>

#define GLM_FORCE_EXPLICIT_CTOR 1
#include <glm/common.hpp>
#include <glm/geometric.hpp>
#include <glm/gtc/matrix_access.hpp>

////////////////////////////////////////////////////////////////////////////////
// Select which standard library
////////////////////////////////////////////////////////////////////////////////

#ifdef BHC_BUILD_CUDA
#include "cuda_runtime.h"
#include "util/UtilsCUDA.cuh"
#define HOST_DEVICE __host__ __device__
// Requires the version of libcudacxx in CUDA 11.5 or later; use the newest CUDA
// version available

// Fixes for some libcudacxx issues for MSVC; won't be needed once a new enough
// version is included in CUDA releases
#if defined(_MSC_VER)

#if _MSC_VER >= 1930
#include <../crt/src/stl/xmath.hpp>
#endif

// Yes, I'm aware this is not GCC--this is a hack for the buggy libcudacxx MSVC
// support
#define __GCC_ATOMIC_BOOL_LOCK_FREE 2
#define __GCC_ATOMIC_CHAR_LOCK_FREE 2
#define __GCC_ATOMIC_CHAR16_T_LOCK_FREE 2
#define __GCC_ATOMIC_CHAR32_T_LOCK_FREE 2
#define __GCC_ATOMIC_WCHAR_T_LOCK_FREE 2
#define __GCC_ATOMIC_SHORT_LOCK_FREE 2
#define __GCC_ATOMIC_INT_LOCK_FREE 2
#define __GCC_ATOMIC_LONG_LOCK_FREE 2
#define __GCC_ATOMIC_LLONG_LOCK_FREE 2
#define __GCC_ATOMIC_POINTER_LOCK_FREE 2

#ifndef __ATOMIC_RELAXED
#define __ATOMIC_RELAXED 0
#define __ATOMIC_CONSUME 1
#define __ATOMIC_ACQUIRE 2
#define __ATOMIC_RELEASE 3
#define __ATOMIC_ACQ_REL 4
#define __ATOMIC_SEQ_CST 5
#endif //__ATOMIC_RELAXED

#endif

#include <cuda/std/complex>
#include <cuda/std/cfloat>
#include <cuda/std/atomic>
#define STD cuda::std
#define BHC_PROGRAMNAME "bellhopcuda"
#else
#define HOST_DEVICE
#include <atomic>
#define STD std
#define BHC_PROGRAMNAME "bellhopcxx"
#endif

////////////////////////////////////////////////////////////////////////////////
// External headers
////////////////////////////////////////////////////////////////////////////////

#include <bhc/bhc.hpp>

namespace bhc {

////////////////////////////////////////////////////////////////////////////////
// Assertions and debug
////////////////////////////////////////////////////////////////////////////////

#define NULLSTATEMENT ((void)0)
#define REQUIRESEMICOLON \
    do { \
        NULLSTATEMENT; \
    } while(false)

// bail() to be used for debugging only, not normal error reporting.
#ifdef __CUDA_ARCH__
[[noreturn]] __forceinline__ __device__ void dev_bail()
{
    __trap();
    __builtin_unreachable();
}
#define bail() dev_bail()
#else
#define bail() throw std::runtime_error("bhc::bail()")
#endif

#ifdef _MSC_VER
#define CPU_NOINLINE __declspec(noinline)
#else
#define CPU_NOINLINE __attribute__((noinline))
#endif

////////////////////////////////////////////////////////////////////////////////
// Real types
////////////////////////////////////////////////////////////////////////////////

#ifdef BHC_USE_FLOATS
#define REAL_MAX FLT_MAX
#define REAL_EPSILON FLT_EPSILON
#define REAL_MINPOS FLT_MIN
#define REAL_PI ((float)M_PI)
#define REAL_REL_SNAP (1.0e-5f)
#else
#define REAL_MAX DBL_MAX
#define REAL_EPSILON DBL_EPSILON
#define REAL_MINPOS DBL_MIN
#define REAL_PI M_PI
#define REAL_REL_SNAP (1.0e-6)
#endif

// BELLHOP uses mostly normal REAL literals, which are float despite most of the
// program using REAL*8 (double). It occasionally uses double-precision literals
// (look like 1.6D-9). For values exactly expressible in both types, e.g. 0.0,
// 1.0, 2.0, 37.0, 0.375, it doesn't matter which type the literal is--except
// that when running in float mode, double-precision literals may cause the
// entire expression to be promoted to double, causing extremely slow double-
// precision instructions to be emitted on CUDA. However, for values not
// expressable as float, e.g. 0.1, the literal type changes the result:
// double d = bar(); float f = foo();
// assert(d * 0.1 == d * 0.1f); //will fail for all d except 0, inf, etc.
// assert((float)(f * 0.1) == f * 0.1f); //will fail for some f
//
// "Real literal"
#ifdef BHC_USE_FLOATS
#define RL(a) (a##f)
#else
#define RL(a) a
#endif
// "Float literal"--This macro is not actually needed, values could just be
// always written as e.g. 0.1f, but this way, every floating-point literal
// should have one or the other macro on it, making it easier to spot errors.
#define FL(a) (a##f)

#define CHECK_REAL_T() \
    static_assert(std::is_floating_point<REAL>::value, "Invalid type for REAL!")

template<typename REAL> REAL fmsub(REAL x, REAL y, REAL z)
{
    // LP: Fused multiply-subtract is a separate instruction in x86. Hopefully
    // this will be emitted, but even if it isn't, the negative z should never
    // change any precision results.
    return STD::fma(x, y, -z);
}

////////////////////////////////////////////////////////////////////////////////
// Complex types
////////////////////////////////////////////////////////////////////////////////

#define J cpx(RL(0.0), RL(1.0))
HOST_DEVICE constexpr inline cpxf Cpx2Cpxf(const cpx &c)
{
    return cpxf((float)c.real(), (float)c.imag());
}
HOST_DEVICE constexpr inline cpx Cpxf2Cpx(const cpxf &c)
{
    return cpx((real)c.real(), (real)c.imag());
}

// CUDA::std::cpx<double> and glm::mat2x2 do not like operators being applied
// with float literals, due to template type deduction issues.
#ifndef BHC_USE_FLOATS
HOST_DEVICE constexpr inline cpx operator-(float a, const cpx &b)
{
    return cpx(a - b.real(), -b.imag());
}
HOST_DEVICE constexpr inline cpx operator*(const cpx &a, float b)
{
    return cpx(a.real() * b, a.imag() * b);
}
HOST_DEVICE constexpr inline cpx operator*(float a, const cpx &b)
{
    return cpx(a * b.real(), a * b.imag());
}
HOST_DEVICE constexpr inline cpx operator/(const cpx &a, float b)
{
    return cpx(a.real() / b, a.imag() / b);
}
HOST_DEVICE inline cpx operator/(float a, const cpx &b) { return (double)a / b; }
HOST_DEVICE inline mat2x2 operator*(float a, const mat2x2 &b) { return (double)a * b; }
#endif

////////////////////////////////////////////////////////////////////////////////
// Misc math
////////////////////////////////////////////////////////////////////////////////

#define SQ(a) ((a) * (a)) // Square
#define CUBE(a) ((a) * (a) * (a))
constexpr real RadDeg = RL(180.0) / REAL_PI;
constexpr real DegRad = REAL_PI / RL(180.0);

template<typename REAL> HOST_DEVICE inline REAL spacing(REAL v)
{
    CHECK_REAL_T();
    return STD::abs(STD::nextafter(v, (REAL)(0.0f)) - v);
}

inline HOST_DEVICE bool isfinite(const vec3 &v)
{
    return STD::isfinite(v.x) && STD::isfinite(v.y) && STD::isfinite(v.z);
}

template<bool X3D> HOST_DEVICE inline constexpr int ZDIM()
{
    if constexpr(X3D)
        return 2;
    else
        return 1;
}
template<typename VEC> HOST_DEVICE inline real &DEP(VEC &v);
template<> HOST_DEVICE inline real &DEP(vec2 &v) { return v.y; }
template<> HOST_DEVICE inline real &DEP(vec3 &v) { return v.z; }
template<typename VEC> HOST_DEVICE inline const real &DEP(const VEC &v);
template<> HOST_DEVICE inline const real &DEP(const vec2 &v) { return v.y; }
template<> HOST_DEVICE inline const real &DEP(const vec3 &v) { return v.z; }
inline HOST_DEVICE vec2 XYCOMP(const vec3 &v) { return vec2(v.x, v.y); }
inline HOST_DEVICE void SETXY(vec3 &v, const vec2 &xy)
{
    v.x = xy.x;
    v.y = xy.y;
}

// max/min are not handled the same way as other math functions by the C++
// standard library and therefore also by libcu++. These versions make sure
// to use the underlying function with the correct precision.
#ifdef __CUDA_ARCH__
#define DEFINE_MATH_FUNC_2(BASE, DEV_F, HOST_F, DEV_D, HOST_D) \
    __device__ inline float BASE(const float &a, const float &b) { return DEV_F(a, b); } \
    __device__ inline double BASE(const double &a, const double &b) \
    { \
        return DEV_D(a, b); \
    }
#define DEFINE_MATH_FUNC_INT_2(BASE, DEV, HOST) \
    __device__ inline int32_t BASE(const int32_t &a, const int32_t &b) \
    { \
        return DEV(a, b); \
    } \
    __device__ inline uint32_t BASE(const uint32_t &a, const uint32_t &b) \
    { \
        return DEV(a, b); \
    } \
    __device__ inline size_t BASE(const size_t &a, const size_t &b) { return DEV(a, b); }
#else
#define DEFINE_MATH_FUNC_2(BASE, DEV_F, HOST_F, DEV_D, HOST_D) \
    inline float BASE(const float &a, const float &b) { return HOST_F(a, b); } \
    inline double BASE(const double &a, const double &b) { return HOST_D(a, b); }
#define DEFINE_MATH_FUNC_INT_2(BASE, DEV, HOST) \
    inline int32_t BASE(const int32_t &a, const int32_t &b) { return HOST(a, b); } \
    inline uint32_t BASE(const uint32_t &a, const uint32_t &b) { return HOST(a, b); } \
    inline size_t BASE(const size_t &a, const size_t &b) { return HOST(a, b); }
#endif

DEFINE_MATH_FUNC_2(max, fmaxf, std::max, fmax, std::max)
DEFINE_MATH_FUNC_2(min, fminf, std::min, fmin, std::min)
DEFINE_MATH_FUNC_INT_2(max, ::max, std::max)
DEFINE_MATH_FUNC_INT_2(min, ::min, std::min)

HOST_DEVICE inline float RealBitsAddInt(float r, int32_t i)
{
#ifdef __CUDA_ARCH__
    return __int_as_float(__float_as_int(r) + i);
#else
    int32_t k;
    std::memcpy(&k, &r, 4);
    k += i;
    float x;
    std::memcpy(&x, &k, 4);
    return x;
#endif
}
HOST_DEVICE inline double RealBitsAddInt(double r, int32_t i)
{
#ifdef __CUDA_ARCH__
    return __longlong_as_double(__double_as_longlong(r) + (int64_t)i);
#else
    int64_t k;
    std::memcpy(&k, &r, 8);
    k += i;
    double x;
    std::memcpy(&x, &k, 8);
    return x;
#endif
}

////////////////////////////////////////////////////////////////////////////////
// Indexing
////////////////////////////////////////////////////////////////////////////////

template<bool O3D> HOST_DEVICE inline int32_t GetNumJobs(
    const Position *Pos, const AnglesStructure *Angles)
{
    int32_t ret = 1;
    if(Angles->alpha.iSingle == 0) ret *= Angles->alpha.n;
    if constexpr(O3D) {
        if(Angles->beta.iSingle == 0) ret *= Angles->beta.n;
        ret *= Pos->NSy;
        ret *= Pos->NSx;
    }
    ret *= Pos->NSz;
    return ret;
}

/**
 * Returns whether the job should continue.
 * `is` changed to `isrc` because `is` is used for steps
 */
template<bool O3D> HOST_DEVICE inline bool GetJobIndices(
    RayInitInfo &rinit, int32_t job, const Position *Pos, const AnglesStructure *Angles)
{
    if(Angles->alpha.iSingle >= 1) {
        // iSingle is 1-indexed because how defined in env file
        rinit.ialpha = Angles->alpha.iSingle - 1;
    } else {
        rinit.ialpha = job % Angles->alpha.n;
        job /= Angles->alpha.n;
    }
    if constexpr(O3D) {
        if(Angles->beta.iSingle >= 1) {
            rinit.ibeta = Angles->beta.iSingle - 1;
        } else {
            rinit.ibeta = job % Angles->beta.n;
            job /= Angles->beta.n;
        }
        rinit.isy = job % Pos->NSy;
        job /= Pos->NSy;
        rinit.isx = job % Pos->NSx;
        job /= Pos->NSx;
    } else {
        rinit.isx = rinit.isy = rinit.ibeta = 0;
    }
    rinit.isz = job;
    return (rinit.isz < Pos->NSz);
}

HOST_DEVICE inline size_t GetFieldAddr(
    int32_t isx, int32_t isy, int32_t isz, int32_t itheta, int32_t id, int32_t ir,
    const Position *Pos)
{
    // clang-format off
    return (((((size_t)isz
        * (size_t)Pos->NSx + (size_t)isx)
        * (size_t)Pos->NSy + (size_t)isy)
        * (size_t)Pos->Ntheta + (size_t)itheta)
        * (size_t)Pos->NRz_per_range + (size_t)id)
        * (size_t)Pos->NRr + (size_t)ir;
    // clang-format on
}

} // namespace bhc

#define _BHC_INCLUDING_COMPONENTS_ 1
#include "util/errors.hpp"
#include "util/prtfileemu.hpp"
#include "util/timing.hpp"
#include "runtype.hpp"
#undef _BHC_INCLUDING_COMPONENTS_

namespace bhc {

////////////////////////////////////////////////////////////////////////////////
// Internal
////////////////////////////////////////////////////////////////////////////////

struct bhcInternal {
    void (*outputCallback)(const char *message);
    std::string FileRoot;
    PrintFileEmu PRTFile;
    std::atomic<int32_t> sharedJobID;
    int gpuIndex, d_multiprocs; // d_warp, d_maxthreads
    int32_t numThreads;
    size_t maxMemory;
    size_t usedMemory;
    bool useRayCopyMode;
    bool noEnvFil;
    uint8_t dim;

    bhcInternal(const bhcInit &init, bool o3d, bool r3d)
        : outputCallback(init.outputCallback),
          FileRoot(
              init.FileRoot == nullptr ? "error_incorrect_use_of_" BHC_PROGRAMNAME
                                       : init.FileRoot),
          PRTFile(this, this->FileRoot, init.prtCallback), gpuIndex(init.gpuIndex),
          numThreads(ModifyNumThreads(init.numThreads)), maxMemory(init.maxMemory),
          usedMemory(0), useRayCopyMode(init.useRayCopyMode),
          noEnvFil(init.FileRoot == nullptr), dim(r3d       ? 3
                                                      : o3d ? 4
                                                            : 2)
    {}
};

template<bool O3D> inline bhcInternal *GetInternal(const bhcParams<O3D> &params)
{
    return reinterpret_cast<bhcInternal *>(params.internal);
}

} // namespace bhc
