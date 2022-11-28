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

////////////////////////////////////////////////////////////////////////////////
// Build configuration
////////////////////////////////////////////////////////////////////////////////

// BHC_DIMMODE is undefined or 0 for "all enabled, set on command line"
#define BHC_ENABLE_2D (!(BHC_DIMMODE == 3 || BHC_DIMMODE == 4))
#define BHC_ENABLE_3D (!(BHC_DIMMODE == 2 || BHC_DIMMODE == 4))
#define BHC_ENABLE_NX2D (!(BHC_DIMMODE == 2 || BHC_DIMMODE == 3))
#if BHC_DIMMODE == 2
#define BHC_DIMNAME "2d"
#elif BHC_DIMMODE == 3
#define BHC_DIMNAME "3d"
#elif BHC_DIMMODE == 4
#define BHC_DIMNAME "nx2d"
#else
#define BHC_DIMNAME
#endif

////////////////////////////////////////////////////////////////////////////////
// General headers
////////////////////////////////////////////////////////////////////////////////

#define _USE_MATH_DEFINES 1 // must be before anything which includes math.h
#include <math.h>

#ifdef _MSC_VER
#define NOMINMAX
#include <windows.h>
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <iomanip>
#include <cctype>
#include <cstring>
#include <string>
#include <locale>
#include <algorithm>
//#include <cfenv>
#include <chrono>
#include <exception>

#include <glm/common.hpp>
#include <glm/geometric.hpp>
#include <glm/gtc/matrix_access.hpp>

////////////////////////////////////////////////////////////////////////////////
// Select which standard library
////////////////////////////////////////////////////////////////////////////////

#ifdef BHC_BUILD_CUDA
#include "cuda_runtime.h"
#define HOST_DEVICE __host__ __device__
// Requires libcu++ 1.4 or higher, which is included in the normal CUDA Toolkit
// install since somewhere between CUDA 11.3 and 11.5. (I.e. 11.5 and later has
// it, 11.2 and earlier does not, not sure about 11.3 and 11.4. You can check by
// seeing if <cuda_install_dir>/include/cuda/std/complex exists.)
#include <cuda/std/complex>
#include <cuda/std/cfloat>
#define STD cuda::std
#define BHC_PROGRAMNAME "bellhopcuda" BHC_DIMNAME
#else
#define HOST_DEVICE
#define STD std
#define BHC_PROGRAMNAME "bellhopcxx" BHC_DIMNAME
#endif

////////////////////////////////////////////////////////////////////////////////
// External headers
////////////////////////////////////////////////////////////////////////////////

#include <bhc/bhc.hpp>

namespace bhc {

////////////////////////////////////////////////////////////////////////////////
// Assertions and debug
////////////////////////////////////////////////////////////////////////////////

extern bool api_okay;

#define NULLSTATEMENT ((void)0)
#define REQUIRESEMICOLON \
    do { \
        NULLSTATEMENT; \
    } while(false)
#define IGNORE_UNUSED(x) \
    do { \
        (void)x; \
    } while(false)

/**
 * Returns a pointer to only the last portion of the source filename.
 */
inline const char *SOURCE_FILENAME(const char *file)
{
    static const char *const tag = "/bellhopcuda/";
    static const int taglen      = 13;
    const char *x                = file;
    for(; *x; ++x) {
        int i = 0;
        for(; i < taglen && x[i]; ++i) {
            if(x[i] != tag[i]) break;
        }
        if(i == taglen) break;
    }
    if(*x) return x + taglen;
    return file;
}

#define BASSERT_STR(x) #x
#define BASSERT_XSTR(x) BASSERT_STR(x)
#ifdef __CUDA_ARCH__
[[noreturn]] __forceinline__ __device__ void dev_bail()
{
    __trap();
    __builtin_unreachable();
}
#define bail() dev_bail()
#define BASSERT(statement) \
    if(__builtin_expect(!(statement), 0)) { \
        GlobalLog("Assertion " #statement " failed line " BASSERT_XSTR(__LINE__) "!\n"); \
        __trap(); \
    } \
    REQUIRESEMICOLON
#else
#define bail() throw std::runtime_error("bhc::bail()")
#define BASSERT(statement) \
    if(!(statement)) { \
        throw std::runtime_error( \
            "Assertion \"" #statement "\" failed in " \
            + std::string(SOURCE_FILENAME(__FILE__)) \
            + " line " BASSERT_XSTR(__LINE__) "\n"); \
    } \
    REQUIRESEMICOLON
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
// Must be below abs(bit_cast<float>(0xFEFEFEFEu) == -1.69e38f)
#define DEBUG_LARGEVAL (1.0e30f)
// #define DEBUG_LARGEVAL (1.0e15f)
#else
#define REAL_MAX DBL_MAX
#define REAL_EPSILON DBL_EPSILON
#define REAL_MINPOS DBL_MIN
#define REAL_PI M_PI
#define REAL_REL_SNAP (1.0e-6)
// Must be below abs(bit_cast<double>(0xFEFEFEFEFEFEFEFEull) == -5.31e303)
#define DEBUG_LARGEVAL (1.0e250)
//#define DEBUG_LARGEVAL (1.0e15)
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

template<bool O3D, bool R3D> HOST_DEVICE VEC23<O3D> RayToOceanX(
    const VEC23<R3D> &x, const Origin<O3D, R3D> &org)
{
    static_assert(O3D || !R3D, "2D ocean but 3D rays not allowed!");
    if constexpr(O3D && !R3D) {
        return vec3(org.xs.x + x.x * org.tradial.x, org.xs.y + x.x * org.tradial.y, x.y);
    } else {
        return x;
    }
}
template<bool O3D, bool R3D> HOST_DEVICE VEC23<O3D> RayToOceanT(
    const VEC23<R3D> &t, const Origin<O3D, R3D> &org)
{
    static_assert(O3D || !R3D, "2D ocean but 3D rays not allowed!");
    if constexpr(O3D && !R3D) {
        return vec3(t.x * org.tradial.x, t.x * org.tradial.y, t.y);
    } else {
        return t;
    }
}
template<bool O3D, bool R3D> HOST_DEVICE VEC23<R3D> OceanToRayX(
    const VEC23<O3D> &x, const Origin<O3D, R3D> &org, const VEC23<R3D> &t)
{
    static_assert(O3D || !R3D, "2D ocean but 3D rays not allowed!");
    if constexpr(O3D && !R3D) {
        // LP: Going back and forth through the coordinate transform won't
        // always keep the precise value, so we may have to finesse the floats.
        vec2 x_orig;
        x_orig.y = x.z;
        if(STD::abs(org.tradial.x) >= STD::abs(org.tradial.y)) {
            x_orig.x = (x.x - org.xs.x) / org.tradial.x;
        } else {
            x_orig.x = (x.x - org.xs.y) / org.tradial.y;
        }
        vec3 x_res = RayToOceanX(x_orig, org);
        if(x_res.x == x.x && x_res.y == x.y) {
            // Got lucky--it went through and came back with the same values.
            return x_orig;
        }
        // Try adding or subtracting one ulp.
        vec2 x_try  = x_orig;
        x_try.x     = RealBitsAddInt(x_orig.x, 1);
        vec2 x_res2 = RayToOceanX(x_try, org);
        if(x_res2.x == x.x && x_res2.y == x.y) return x_try;
        x_try.x = RealBitsAddInt(x_orig.x, -1);
        x_res2  = RayToOceanX(x_try, org);
        if(x_res2.x == x.x && x_res2.y == x.y) return x_try;
        // No hope of being exact. Just try to be slightly forward of the boundary.
        x_try.x = x_orig.x + RL(1e-6) * t.x;
        if(x_try.x == x_orig.x) { x_try.x = x_orig.x + RL(1e-3) * t.x; }
        return x_try;
    } else {
        return x;
    }
}

////////////////////////////////////////////////////////////////////////////////
// String manipulation
////////////////////////////////////////////////////////////////////////////////

inline bool isInt(std::string str, bool allowNegative = true)
{
    if(str.empty()) return false;
    for(size_t i = 0; i < str.length(); ++i) {
        if(str[i] == '-') {
            if(i != 0 || !allowNegative || str.length() == 1) return false;
            continue;
        } else if(str[i] >= '0' && str[i] <= '9') {
            continue;
        } else {
            return false;
        }
    }
    return true;
}

inline bool isReal(std::string str)
{
    if(str.empty()) return false;
    char *ptr;
    strtod(str.c_str(), &ptr);
    return (*ptr) == '\0';
}

inline std::ostream &operator<<(std::ostream &s, const cpx &v)
{
    s << "(" << v.real() << "," << v.imag() << ")";
    return s;
}

inline std::ostream &operator<<(std::ostream &s, const vec2 &v)
{
    s << v.x << " " << v.y;
    return s;
}

// Courtesy Evan Teran, https://stackoverflow.com/a/217605
// trim from start (in place)
static inline void ltrim(std::string &s)
{
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
                return !std::isspace(ch);
            }));
}
// trim from end (in place)
static inline void rtrim(std::string &s)
{
    s.erase(
        std::find_if(
            s.rbegin(), s.rend(), [](unsigned char ch) { return !std::isspace(ch); })
            .base(),
        s.end());
}
// trim from both ends (in place)
static inline void trim(std::string &s)
{
    ltrim(s);
    rtrim(s);
}
// trim from start (copying)
static inline std::string ltrim_copy(std::string s)
{
    ltrim(s);
    return s;
}
// trim from end (copying)
static inline std::string rtrim_copy(std::string s)
{
    rtrim(s);
    return s;
}
// trim from both ends (copying)
static inline std::string trim_copy(std::string s)
{
    trim(s);
    return s;
}

////////////////////////////////////////////////////////////////////////////////
// Other components
////////////////////////////////////////////////////////////////////////////////

} // namespace bhc

#define _BHC_INCLUDING_COMPONENTS_ 1
#include "logging.hpp"
#include "ldio.hpp"
#include "bino.hpp"
#include "prtfileemu.hpp"
#include "atomics.hpp"

#ifdef BHC_BUILD_CUDA
#include "UtilsCUDA.cuh"
#endif

#undef _BHC_INCLUDING_COMPONENTS_

namespace bhc {

////////////////////////////////////////////////////////////////////////////////
// CUDA memory
////////////////////////////////////////////////////////////////////////////////

template<typename T> inline T *allocate(size_t n = 1)
{
    T *ret;
#ifdef BHC_BUILD_CUDA
    checkCudaErrors(cudaMallocManaged(&ret, n * sizeof(T)));
#else
    ret = new T[n];
#endif
#ifdef BHC_DEBUG
    // Debugging: Fill memory with garbage data to help detect uninitialized vars
    memset(ret, 0xFE, n * sizeof(T));
#endif
    return ret;
}

template<typename T> inline void deallocate(T *&ptr)
{
#ifdef BHC_BUILD_CUDA
    checkCudaErrors(cudaFree(ptr));
#else
    delete[] ptr;
#endif
    ptr = nullptr;
}

template<typename T> inline void checkdeallocate(T *&ptr)
{
    if(ptr != nullptr) deallocate(ptr);
}

template<typename T> inline void checkallocate(T *&ptr, size_t n = 1)
{
    if(ptr != nullptr) deallocate(ptr);
    ptr = allocate<T>(n);
}

////////////////////////////////////////////////////////////////////////////////
// Algorithms
////////////////////////////////////////////////////////////////////////////////

template<typename T> inline void Sort(T *arr, size_t n) { std::sort(arr, arr + n); }
template<> inline void Sort(cpx *arr, size_t n)
{
    // Workaround because both libstdc++ and libcu++ provide swap
    std::complex<real> *arr2 = reinterpret_cast<std::complex<real> *>(arr);
    std::sort(
        arr2, arr2 + n, [](const std::complex<real> &a, const std::complex<real> &b) {
            return a.real() > b.real(); // Based on order of decreasing real part
        });
}

/**
 * mbp: tests whether an input vector is strictly monotonically increasing
 */
template<typename REAL> HOST_DEVICE inline bool monotonic(REAL *arr, size_t n)
{
    CHECK_REAL_T();
    if(n == 1) return true;
    for(size_t i = 0; i < n - 1; ++i) {
        if(arr[i + 1] <= arr[i]) return false;
    }
    return true;
}
template<typename VEC2> HOST_DEVICE inline bool monotonic(
    VEC2 *arr, size_t n, const int32_t dim)
{
    if(n == 1) return true;
    for(size_t i = 0; i < n - 1; ++i) {
        if(arr[i + 1][dim] <= arr[i][dim]) return false;
    }
    return true;
}
template<typename REAL> HOST_DEVICE inline bool monotonic(
    REAL *arr, size_t n, const int32_t stridereals, const int32_t offset)
{
    CHECK_REAL_T();
    if(n == 1) return true;
    for(size_t i = 0; i < n - 1; ++i) {
        if(arr[(i + 1) * stridereals + offset] <= arr[i * stridereals + offset])
            return false;
    }
    return true;
}

/**
 * Returns the index of the largest element in the array less than or equal to
 * target, given that the array is monotonically increasing.
 * It is not a bug that arr is REAL and target is real. If the array is of
 * reduced precision relative to the actual desired target, the comparisons
 * must be done in the higher precision, or incorrect results can be returned
 * (e.g. round(target) <= arr[i] but target > arr[i]).
 */
template<typename REAL> HOST_DEVICE inline int32_t BinarySearchLEQ(
    REAL *arr, int32_t n, const int32_t stride, const int32_t offset, real target)
{
    CHECK_REAL_T();
    int32_t low = 0; // Low is included
    int32_t hi  = n; // Hi is excluded
    while(low < hi) {
        int32_t t = (low + hi) / 2;
        if((real)arr[t * stride + offset] > target) {
            hi = t;
        } else if(t >= n - 1) {
            return n - 1;
        } else if((real)arr[(t + 1) * stride + offset] > target) {
            return t;
        } else {
            low = t + 1;
        }
    }
    return low;
}

/**
 * Returns the index of the smallest element in the array greater than or equal
 * to target, given that the array is monotonically increasing.
 * It is not a bug that arr is REAL and target is real. If the array is of
 * reduced precision relative to the actual desired target, the comparisons
 * must be done in the higher precision, or incorrect results can be returned
 * (e.g. round(target) >= arr[i] but target < arr[i]).
 */
template<typename REAL> HOST_DEVICE inline int32_t BinarySearchGEQ(
    REAL *arr, int32_t n, const int32_t stride, const int32_t offset, real target)
{
    CHECK_REAL_T();
    int32_t low = -1;    // Low is excluded
    int32_t hi  = n - 1; // Hi is included
    while(low < hi) {
        int32_t t = (low + hi + 1) / 2; // Round up
        if((real)arr[t * stride + offset] < target) {
            low = t;
        } else if(t <= 0) {
            return 0;
        } else if((real)arr[(t - 1) * stride + offset] < target) {
            return t;
        } else {
            hi = t - 1;
        }
    }
    return hi;
}

/**
 * Returns the index of the smallest element in the array strictly greater than
 * to target, given that the array is monotonically increasing.
 * It is not a bug that arr is REAL and target is real. If the array is of
 * reduced precision relative to the actual desired target, the comparisons
 * must be done in the higher precision, or incorrect results can be returned
 * (e.g. round(target) > arr[i] but target <= arr[i]).
 */
template<typename REAL> HOST_DEVICE inline int32_t BinarySearchGT(
    REAL *arr, int32_t n, const int32_t stride, const int32_t offset, real target)
{
    CHECK_REAL_T();
    int32_t low = -1;    // Low is excluded
    int32_t hi  = n - 1; // Hi is included
    while(low < hi) {
        int32_t t = (low + hi + 1) / 2; // Round up
        if((real)arr[t * stride + offset] <= target) {
            low = t;
        } else if(t <= 0) {
            return 0;
        } else if((real)arr[(t - 1) * stride + offset] <= target) {
            return t;
        } else {
            hi = t - 1;
        }
    }
    return hi;
}

/**
 * mbp: full 360-degree sweep? remove duplicate angle/beam
 */
template<typename REAL> HOST_DEVICE inline void CheckFix360Sweep(
    const REAL *angles, int32_t &n)
{
    CHECK_REAL_T();
    // LP: Changed from (the FORTRAN equivalent of) REAL_MINPOS, see Fortran
    // version readme.
    if(n > 1
       && STD::abs(STD::fmod(angles[n - 1] - angles[0], FL(360.0)))
           < FL(10.0) * spacing(RL(360.0)))
        --n;
}

template<typename REAL> inline void EchoVector(
    REAL *v, int32_t Nv, PrintFileEmu &PRTFile, int32_t NEcho = 10,
    const char *ExtraSpaces = "")
{
    CHECK_REAL_T();
    PRTFile << std::setprecision(6) << ExtraSpaces;
    for(int32_t i = 0, r = 0; i < bhc::min(Nv, NEcho); ++i) {
        PRTFile << std::setw(14) << v[i] << " ";
        ++r;
        if(r == 5) {
            r = 0;
            PRTFile << "\n" << ExtraSpaces;
        }
    }
    if(Nv > NEcho) PRTFile << "... " << std::setw(14) << v[Nv - 1];
    PRTFile << "\n";
    /*
    PRTFile << std::setprecision(12);
    for(int32_t i=0; i<Nv; ++i) PRTFile << std::setw(20) << v[i] << "\n";
    PRTFile << "\n";
    */
}

/**
 * If x[2] == -999.9 then subtabulation is performed
 * i.e., a vector is generated with Nx points in [x[0], x[1]]
 * If x[1] == -999.9 then x[0] is repeated into x[1]
 */
template<typename REAL> HOST_DEVICE inline void SubTab(REAL *x, int32_t Nx)
{
    CHECK_REAL_T();
    // LP: Must be literals in the type of REAL. (double)(0.01f) != 0.01 and
    // (float)(0.01) is not guaranteed to equal 0.01f.
    // (For proof, we write the true real number as 0.01* and double as 0.01d.
    // Suppose 0.01f = 0.01* + 4.9e-9* because it is the closest floating point
    // value to 0.01*, because the next number below is 0.01* - 5.1e-9.
    // Now suppose 0.01d = 0.01* - 0.2e-9 (because the next number up is
    // 0.01* + 0.3e-9 or whatever). Now, (float)(0.01d) will evaluate to
    // 0.01* - 5.1e-9, because 0.01* - 0.2e-9 is closer to 0.01* - 5.1e-9 than
    // it is to 0.01* + 4.9e-9*.)
    REAL minus999, pointohone;
    if constexpr(sizeof(REAL) == 4) {
        minus999   = -999.9f;
        pointohone = 0.01f;
    } else {
        minus999   = -999.9;
        pointohone = 0.01;
    }
    if(Nx >= 3) {
        if(STD::abs(x[2] - minus999) < pointohone) {
            if(STD::abs(x[1] - minus999) < pointohone) x[1] = x[0];
            REAL deltax = (x[1] - x[0]) / (REAL)(Nx - 1);
            REAL x0     = x[0];
            for(int32_t i = 0; i < Nx; ++i) {
                // LP: gfortran emits AVX-512 FMA instructions here. Without an
                // FMA, the resulting receiver range positions are slightly
                // different at bondaries, leading to inconsistencies later
                // x[i] = x0 + (REAL)i * deltax;
                x[i] = STD::fma((REAL)i, deltax, x0);
            }
        }
    }
}

HOST_DEVICE inline void RayNormalImpl(
    const vec3 &t, real phi, bool ignorephi0, real c, vec3 &e1, vec3 &e2)
{
    real rl = glm::length(vec2(t.x, t.y));

    if(phi != RL(0.0) || ignorephi0) {
        real cosphi = STD::cos(phi), sinphi = STD::sin(phi);

        // e1
        e1.x = (c * t.x * t.z * cosphi + t.y * sinphi) / rl;
        e1.y = (c * t.y * t.z * cosphi - t.x * sinphi) / rl;
        e1.z = -c * rl * cosphi;

        // e2
        e2.x = (c * t.x * t.z * sinphi - t.y * cosphi) / rl;
        e2.y = (c * t.y * t.z * sinphi + t.x * cosphi) / rl;
        e2.z = -c * rl * sinphi;
        /*
        // LP: This algorithm is sensitive to precision. e1 is used to compute
        // an update to phi, and phi is used here to compute e1. Even if the
        // other variables (t, c, rayn1) do not diverge between Fortran and C++,
        // phi and e1 may diverge due to floating-point precision.
        // This implementation matches the set of AVX-512 operations performed
        // by the gfortran build of BELLHOP.
        real t1t3 = t.x * t.z;
        real t2t3 = t.y * t.z;
        real ct1t3 = c * t1t3;
        real ct2t3 = c * t2t3;
        real crl   = c * rl;
        real rcprl = RL(1.0) / rl;
        real e11build = sinphi * t.y;
        real e12build = sinphi * t.x;
        e11build = STD::fma(cosphi, ct1t3, e11build);
        real t2cosphi = cosphi * t.y;
        e12build =    fmsub(cosphi, ct2t3, e12build);
        real e13build = cosphi * crl;
        real e22build = cosphi * t.x;
        real e21build = fmsub(ct1t3, sinphi, t2cosphi);
        e1.z = -e13build;
        e22build = STD::fma(sinphi, ct2t3, e22build);
        real e23build = sinphi * crl;
        e1.x = e11build * rcprl;
        e1.y = e12build * rcprl;
        e2.z = -e23build;
        e2.x = e21build * rcprl;
        e2.y = e22build * rcprl;
        */
    } else {
        e1 = vec3(c * t.x * t.z / rl, c * t.y * t.z / rl, -c * rl);
        e2 = vec3(-t.y / rl, t.x / rl, RL(0.0));
    }
}
/**
 * Computes the ray normals
 *
 * t: tangent vector (NOT) normalized
 * phi: torsion
 * c: sound speed
 * e1, e2: ray unit normals
 */
HOST_DEVICE inline void RayNormal(const vec3 &t, real phi, real c, vec3 &e1, vec3 &e2)
{
    RayNormalImpl(t, phi, false, c, e1, e2);
}
/**
 * Computes the ray normals
 * Same as routine RayNormal except this version assumes t is already normalized
 */
HOST_DEVICE inline void RayNormal_unit(const vec3 &t, real phi, vec3 &e1, vec3 &e2)
{
    RayNormalImpl(t, phi, true, RL(1.0), e1, e2);
}

////////////////////////////////////////////////////////////////////////////////
// Timing
////////////////////////////////////////////////////////////////////////////////

class Stopwatch {
public:
    Stopwatch() {}
    inline void tick() { tstart = std::chrono::high_resolution_clock::now(); }
    inline void tock()
    {
        using namespace std::chrono;
        high_resolution_clock::time_point tend = high_resolution_clock::now();
        double dt = (duration_cast<duration<double>>(tend - tstart)).count();
        dt *= 1000.0;
        GlobalLog("%f ms\n", dt);
    }

private:
    std::chrono::high_resolution_clock::time_point tstart;
};

HOST_DEVICE inline void PrintMatrix(const mat2x2 &m, const char *label)
{
    GlobalLog(
        "%s: /%12.7e %12.7e\\\n       \\%12.7e %12.7e/\n", label, m[0][0], m[1][0],
        m[0][1], m[1][1]);
}

} // namespace bhc
