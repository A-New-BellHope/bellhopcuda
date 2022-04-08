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

////////////////////////////////////////////////////////////////////////////////
//Select which standard library
////////////////////////////////////////////////////////////////////////////////

#ifdef BHC_BUILD_CUDA
#include "cuda_runtime.h"
#include "UtilsCUDA.cuh"
#define HOST_DEVICE __host__ __device__
//Requires libcu++ 1.4 or higher, which is included in the normal CUDA Toolkit
//install since somewhere between CUDA 11.3 and 11.5. (I.e. 11.5 and later has
//it, 11.2 and earlier does not, not sure about 11.3 and 11.4. You can check by
//seeing if <cuda_install_dir>/include/cuda/std/complex exists.)
#include <cuda/std/complex>
#include <cuda/std/cfloat>
#define STD cuda::std
#define BHC_PROGRAMNAME "bellhopcuda"
#else
#define HOST_DEVICE
#define STD std
#define BHC_PROGRAMNAME "bellhopcxx"
#endif

////////////////////////////////////////////////////////////////////////////////
//External headers
////////////////////////////////////////////////////////////////////////////////

#include <bhc/bhc.hpp>

namespace bhc {

////////////////////////////////////////////////////////////////////////////////
//Assertions and debug
////////////////////////////////////////////////////////////////////////////////

extern bool api_okay;

#define NULLSTATEMENT ((void)0)
#define REQUIRESEMICOLON do{NULLSTATEMENT;} while(false)
#define IGNORE_UNUSED(x) do{ (void)x; } while(false)

/**
 * Returns a pointer to only the last portion of the source filename.
*/
inline const char *SOURCE_FILENAME(const char *file){
    static const char *const tag = "/bellhopcuda/";
    static const int taglen = 13;
    const char *x = file;
    for(; *x; ++x){
        int i=0;
        for(; i<taglen && x[i]; ++i){
            if(x[i] != tag[i]) break;
        }
        if(i==taglen) break;
    }
    if(*x) return x+taglen;
    return file;
}

#define BASSERT_STR(x) #x
#define BASSERT_XSTR(x) BASSERT_STR(x)
#ifdef __CUDA_ARCH__
#define bail __trap
#define BASSERT(statement) \
if(__builtin_expect(!(statement), 0)) { \
	printf("Assertion " #statement " failed line " BASSERT_XSTR(__LINE__) "!\n"); \
	__trap(); \
} REQUIRESEMICOLON
#else
#define bail() throw std::runtime_error("bhc::bail()")
#define BASSERT(statement) \
if(!(statement)){ \
	throw std::runtime_error("Assertion \"" #statement "\" failed in " \
		+ std::string(SOURCE_FILENAME(__FILE__)) + " line " BASSERT_XSTR(__LINE__) "\n"); \
} REQUIRESEMICOLON
#endif

////////////////////////////////////////////////////////////////////////////////
//Real types
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
//Must be below abs(bit_cast<double>(0xFEFEFEFEFEFEFEFEull) == -5.31e303)
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

////////////////////////////////////////////////////////////////////////////////
//Complex types
////////////////////////////////////////////////////////////////////////////////

#define J cpx(RL(0.0), RL(1.0))
HOST_DEVICE constexpr inline cpxf Cpx2Cpxf(const cpx &c){
	return cpxf((float)c.real(), (float)c.imag());
}
HOST_DEVICE constexpr inline cpx Cpxf2Cpx(const cpxf &c){
	return cpx((real)c.real(), (real)c.imag());
}

//CUDA::std::cpx<double> does not like operators being applied with float
//literals, due to template type deduction issues.
#ifndef BHC_USE_FLOATS
HOST_DEVICE constexpr inline cpx operator-(float a, const cpx &b){
    return cpx(a - b.real(), -b.imag());
}
HOST_DEVICE constexpr inline cpx operator*(const cpx &a, float b){
    return cpx(a.real() * b, a.imag() * b);
}
HOST_DEVICE constexpr inline cpx operator*(float a, const cpx &b){
    return cpx(a * b.real(), a * b.imag());
}
HOST_DEVICE constexpr inline cpx operator/(const cpx &a, float b){
    return cpx(a.real() / b, a.imag() / b);
}
HOST_DEVICE inline cpx operator/(float a, const cpx &b){
    return (double)a / b;
}
#endif

////////////////////////////////////////////////////////////////////////////////
//Misc math
////////////////////////////////////////////////////////////////////////////////

#define SQ(a) ((a) * (a)) //Square
#define CUBE(a) ((a) * (a) * (a))
constexpr real RadDeg = RL(180.0) / REAL_PI;
constexpr real DegRad = REAL_PI / RL(180.0);

template<typename REAL> HOST_DEVICE inline REAL spacing(REAL v){
	return STD::abs(STD::nextafter(v, (REAL)(0.0f)) - v);
}

inline HOST_DEVICE bool isfinite(const vec3 &v){
    return STD::isfinite(v.x) && STD::isfinite(v.y) && STD::isfinite(v.z);
}

//max/min are not handled the same way as other math functions by the C++
//standard library and therefore also by libcu++. These versions make sure
//to use the underlying function with the correct precision.
#ifdef __CUDA_ARCH__
#define DEFINE_MATH_FUNC_2(BASE, DEV_F, HOST_F, DEV_D, HOST_D) \
	__device__ inline float BASE(const float &a, const float &b) { return DEV_F(a, b); } \
	__device__ inline double BASE(const double &a, const double &b) { return DEV_D(a, b); }
#define DEFINE_MATH_FUNC_INT_2(BASE, DEV, HOST) \
	__device__ inline int32_t BASE(const int32_t &a, const int32_t &b) { return DEV(a, b); } \
	__device__ inline uint32_t BASE(const uint32_t &a, const uint32_t &b) { return DEV(a, b); } \
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

////////////////////////////////////////////////////////////////////////////////
//String manipulation
////////////////////////////////////////////////////////////////////////////////

inline bool isInt(std::string str, bool allowNegative = true){
	if(str.empty()) return false;
	for(size_t i=0; i<str.length(); ++i){
		if(str[i] == '-'){
			if(i != 0 || !allowNegative || str.length() == 1) return false;
			continue;
		}else if(str[i] >= '0' && str[i] <= '9'){
			continue;
		}else{
			return false;
		}
	}
	return true;
}

inline bool isReal(std::string str){
	if(str.empty()) return false;
	char *ptr;
	strtod(str.c_str(), &ptr);
	return (*ptr) == '\0';
}

inline std::ostream &operator<<(std::ostream &s, const cpx &v){
	s << "(" << v.real() << "," << v.imag() << ")";
	return s;
}

inline std::ostream &operator<<(std::ostream &s, const vec2 &v){
	s << v.x << " " << v.y;
	return s;
}

//Courtesy Evan Teran, https://stackoverflow.com/a/217605
// trim from start (in place)
static inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}
// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}
// trim from both ends (in place)
static inline void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}
// trim from start (copying)
static inline std::string ltrim_copy(std::string s) {
    ltrim(s);
    return s;
}
// trim from end (copying)
static inline std::string rtrim_copy(std::string s) {
    rtrim(s);
    return s;
}
// trim from both ends (copying)
static inline std::string trim_copy(std::string s) {
    trim(s);
    return s;
}

////////////////////////////////////////////////////////////////////////////////
//Other components
////////////////////////////////////////////////////////////////////////////////

}

#define _BHC_INCLUDING_COMPONENTS_ 1
#include "ldio.hpp"
#include "bino.hpp"
#include "prtfileemu.hpp"
#include "atomics.hpp"
#undef _BHC_INCLUDING_COMPONENTS_

namespace bhc {

////////////////////////////////////////////////////////////////////////////////
//CUDA memory
////////////////////////////////////////////////////////////////////////////////

template<typename T> inline T* allocate(size_t n=1){
	T* ret;
	#ifdef BHC_BUILD_CUDA
		checkCudaErrors(cudaMallocManaged(&ret, n*sizeof(T)));
	#else
		ret = new T[n];
	#endif
	// Debugging: Fill memory with garbage data to help detect uninitialized vars
	memset(ret, 0xFE, n*sizeof(T));
	return ret;
}

template<typename T> inline void deallocate(T *&ptr){
	#ifdef BHC_BUILD_CUDA
		checkCudaErrors(cudaFree(ptr));
	#else
		delete[] ptr;
	#endif
	ptr = nullptr;
}

template<typename T> inline void checkdeallocate(T *&ptr){
    if(ptr != nullptr) deallocate(ptr);
}

template<typename T> inline void checkallocate(T *&ptr, size_t n=1){
    if(ptr != nullptr) deallocate(ptr);
    ptr = allocate<T>(n);
}

////////////////////////////////////////////////////////////////////////////////
//Algorithms
////////////////////////////////////////////////////////////////////////////////

template<typename T> inline void Sort(T *arr, size_t n){
	std::sort(arr, arr + n);
}
template<> inline void Sort(cpx *arr, size_t n){
    //Workaround because both libstdc++ and libcu++ provide swap
    std::complex<real> *arr2 = reinterpret_cast<std::complex<real> *>(arr);
	std::sort(arr2, arr2 + n, [](const std::complex<real> &a, const std::complex<real> &b){
		return a.real() > b.real(); // Based on order of decreasing real part
	});
}

/**
 * mbp: tests whether an input vector is strictly monotonically increasing
 */
template<typename REAL> HOST_DEVICE inline bool monotonic(REAL *arr, size_t n){
	if(n == 1) return true;
	for(size_t i=0; i<n-1; ++i){
		if(arr[i+1] <= arr[i]) return false;
	}
	return true;
}
template<typename VEC2> HOST_DEVICE inline bool monotonic(VEC2 *arr, size_t n,
	const int32_t dim){
	if(n == 1) return true;
	for(size_t i=0; i<n-1; ++i){
		if(arr[i+1][dim] <= arr[i][dim]) return false;
	}
	return true;
}
template<typename REAL> HOST_DEVICE inline bool monotonic(REAL *arr, size_t n,
	const int32_t stridereals, const int32_t offset)
{
	if(n == 1) return true;
	for(size_t i=0; i<n-1; ++i){
		if(arr[(i+1)*stridereals+offset] <= arr[i*stridereals+offset]){
            printf("%g <= %g\n", arr[(i+1)*stridereals+offset], arr[i*stridereals+offset]);
            return false;
        }
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
template<typename REAL> HOST_DEVICE inline int32_t BinarySearchLEQ(REAL *arr, 
	int32_t n, const int32_t stride, const int32_t offset, real target)
{
	int32_t low = 0; //Low is included
	int32_t hi = n; //Hi is excluded
	while(low < hi){
		int32_t t = (low + hi) / 2;
		if((real)arr[t*stride+offset] > target){
			hi = t;
		}else if(t >= n-1){
			return n-1;
		}else if((real)arr[(t+1)*stride+offset] > target){
			return t;
		}else{
			low = t+1;
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
template<typename REAL> HOST_DEVICE inline int32_t BinarySearchGEQ(REAL *arr,
	int32_t n, const int32_t stride, const int32_t offset, real target)
{
	int32_t low = -1; //Low is excluded
	int32_t hi = n-1; //Hi is included
	while(low < hi){
		int32_t t = (low + hi + 1) / 2; //Round up
		if((real)arr[t*stride+offset] < target){
			low = t;
		}else if(t <= 0){
			return 0;
		}else if((real)arr[(t-1)*stride+offset] < target){
			return t;
		}else{
			hi = t-1;
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
template<typename REAL> HOST_DEVICE inline int32_t BinarySearchGT(REAL *arr,
	int32_t n, const int32_t stride, const int32_t offset, real target)
{
	int32_t low = -1; //Low is excluded
	int32_t hi = n-1; //Hi is included
	while(low < hi){
		int32_t t = (low + hi + 1) / 2; //Round up
		if((real)arr[t*stride+offset] <= target){
			low = t;
		}else if(t <= 0){
			return 0;
		}else if((real)arr[(t-1)*stride+offset] <= target){
			return t;
		}else{
			hi = t-1;
		}
	}
	return hi;
}

/**
 * mbp: full 360-degree sweep? remove duplicate angle/beam
 */
template<typename REAL> HOST_DEVICE inline void CheckFix360Sweep(const REAL *angles,
	int32_t &n){
	if(n > 1 && STD::abs(STD::fmod(angles[n-1] - angles[0], FL(360.0)))
            < FL(10.0) * spacing(RL(1.0)))
        --n;
}

template<typename REAL> inline void EchoVector(REAL *v, int32_t Nv,
    PrintFileEmu &PRTFile, int32_t NEcho = 10)
{
    PRTFile << std::setprecision(6);
    for(int32_t i=0, r=0; i<bhc::min(Nv, NEcho); ++i){
        PRTFile << std::setw(14) << v[i] << " ";
        ++r;
        if(r == 5){
            r = 0;
            PRTFile << "\n";
        }
    }
    if(Nv > NEcho) PRTFile << "... " << std::setw(14) << v[Nv-1];
    PRTFile << "\n";
    /*
    PRTFile << std::setprecision(12);
    for(int32_t i=0; i<Nv; ++i) PRTFile << std::setw(20) << v[i] << "\n";
    PRTFile << "\n";
    */
}

template<typename REAL> HOST_DEVICE inline void SubTab(REAL *x, int32_t Nx)
{
    if(Nx >= 3){
        if(STD::abs(x[2] - (REAL)(-999.9f)) < (REAL)(0.01f)){ 
            if(STD::abs(x[1] - (REAL)(-999.9f)) < (REAL)(0.01f)) x[1] = x[0];
            REAL deltax = (x[1] - x[0]) / (REAL)(Nx - 1);
            REAL x0 = x[0];
            for(int32_t i=0; i<Nx; ++i){
                // LP: gfortran emits AVX-512 FMA instructions here. Without an
                // FMA, the resulting receiver range positions are slightly
                // different at bondaries, leading to inconsistencies later
                //x[i] = x0 + (REAL)i * deltax;
                x[i] = STD::fma((REAL)i, deltax, x0);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
//Timing
////////////////////////////////////////////////////////////////////////////////

class Stopwatch
{
public:
    Stopwatch() {}
    inline void tick() {
        tstart = std::chrono::high_resolution_clock::now();
    }
    inline void tock() {
        using namespace std::chrono;
        high_resolution_clock::time_point tend = high_resolution_clock::now();
        double dt = (duration_cast<duration<double>>(tend - tstart)).count();
        dt *= 1000.0;
        std::cout << dt << " ms\n";
    }
private:
    std::chrono::high_resolution_clock::time_point tstart;
};

}
