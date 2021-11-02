#pragma once

////////////////////////////////////////////////////////////////////////////////
//Common includes
////////////////////////////////////////////////////////////////////////////////

#define _USE_MATH_DEFINES 1 //must be before anything which includes math.h
#include <math.h>
#include <algorithm>
#include <locale>
#include <cctype>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <string>

////////////////////////////////////////////////////////////////////////////////
//Select which standard library
////////////////////////////////////////////////////////////////////////////////

#ifdef BUILD_CUDA
#define HOST_DEVICE __host__ __device__
#include <cuda/std/complex>
#include <cuda/std/cfloat>
//libcu++
#define STD cuda::std
#else
#define HOST_DEVICE
#include <complex>
#include <cfloat>
#define STD std
#endif

#ifdef __CUDA_ARCH__
#define bail __trap
#else
#define bail std::abort
#endif

////////////////////////////////////////////////////////////////////////////////
//Math types
////////////////////////////////////////////////////////////////////////////////

#ifdef USE_FLOATS
using real = float;
#define REAL_MAX FLT_MAX
#define REAL_EPSILON FLT_EPSILON
//Real constant: transforms 2.0 into 2.0f. This is needed on CUDA or else
//unwanted double-precision instructions may be emitted.
#define RC(a) (a##f)
#else
using real = double;
#define REAL_MAX DBL_MAX
#define REAL_EPSILON DBL_EPSILON
#define RC(a) a
#endif

using cpx = STD::complex<real>;
constexpr cpx J = cpx(RC(0.0), RC(1.0));

#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
//#include <glm/vec4.hpp>
//#include <glm/mat2x2.hpp>
//#include <glm/mat3x3.hpp>
//#include <glm/mat4x4.hpp>
#include <glm/common.hpp>
#include <glm/geometric.hpp>

using vec2 = glm::vec<2, real, glm::defaultp>;
using vec3 = glm::vec<3, real, glm::defaultp>;
// using vec4 = glm::vec<4, real, glm::defaultp>;
// using mat2 = glm::mat<2, 2, real, glm::defaultp>;
// using mat3 = glm::mat<3, 3, real, glm::defaultp>;
// using mat4 = glm::mat<4, 4, real, glm::defaultp>;

#define SQ(a) ((a) * (a)) //Square
#define CUBE(a) ((a) * (a) * (a))
constexpr real RadDeg = RC(180.0) / M_PI;
constexpr real DegRad = M_PI / RC(180.0);

inline bool isInt(std::string str, bool allowNegative = true){
	if(str.empty()) return false;
	for(int i=0; i<str.length(); ++i){
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

////////////////////////////////////////////////////////////////////////////////
//String manipulation
////////////////////////////////////////////////////////////////////////////////

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
//CUDA memory
////////////////////////////////////////////////////////////////////////////////

template<typename T> inline T* allocate(size_t n=1){
	#ifdef BUILD_CUDA
		T* ret;
		checkCudaErrors(cudaMallocManaged(&ret, n*sizeof(T)));
		return ret;
	#else
		return new T[n];
	#endif
}

template<typename T> inline void deallocate(T *&ptr){
	#ifdef BUILD_CUDA
		checkCudaErrors(cudaFree(ptr));
	#else
		delete[] ptr;
	#endif
	ptr = nullptr;
}

////////////////////////////////////////////////////////////////////////////////
//Algorithms
////////////////////////////////////////////////////////////////////////////////

template<typename T> inline void Sort(T *arr, size_t n){
	std::sort(arr, arr + n);
}
template<> inline void Sort(cpx *arr, size_t n){
	std::sort(arr, arr + n, [](const cpx &a, const cpx &b){
		return a.real() > b.real(); // Based on order of decreasing real part
	});
}

/**
 * mbp: tests whether an input vector is strictly monotonically increasing
 */
HOST_DEVICE inline bool monotonic(real *arr, size_t n){
	if(n == 1) return true;
	for(size_t i=0; i<n-1; ++i){
		if(arr[i+1] <= arr[i]) return false;
	}
	return true;
}
HOST_DEVICE inline bool monotonic(vec2 *arr, size_t n, const int32_t dim){
	if(n == 1) return true;
	for(size_t i=0; i<n-1; ++i){
		if(arr[i+1][dim] <= arr[i][dim]) return false;
	}
	return true;
}
HOST_DEVICE inline bool monotonic(real *arr, size_t n, const int32_t stridereals, 
	const int32_t offset)
{
	if(n == 1) return true;
	for(size_t i=0; i<n-1; ++i){
		if(arr[(i+1)*stridereals+offset] <= arr[i*stridereals+offset]) return false;
	}
	return true;
}

/**
 * Returns the index of the element in the array less than or equal to target,
 * given that the array is monotonically increasing.
 */
HOST_DEVICE inline int32_t BinarySearch(real *arr, int32_t n, const int32_t stride,
	const int32_t offset, real target)
{
	int32_t low = 0;
	int32_t hi = n;
	while(low < hi){
		int32_t t = (low + hi) / 2;
		if(arr[t*stride+offset] > target){
			hi = t;
		}else if(t >= n-1){
			return n-1;
		}else if(arr[(t+1)*stride+offset] >= target){
			return t;
		}else{
			low = t+1;
		}
	}
	return low;
}

/**
 * mbp: full 360-degree sweep? remove duplicate angle/beam
 */
HOST_DEVICE inline void CheckFix360Sweep(const real *angles, int32_t &n){
	if(n > 1 && STD::abs(STD::fmod(angles[n-1] - angles[0], RC(360.0)))
            < RC(10.0) * (STD::nextafter(RC(1.0), RC(2.0)) - RC(1.0)))
        --n;
}

inline void EchoVector(real *v, int32_t Nv, std::ofstream &PRTFile)
{
    constexpr int32_t NEcho = 10;
    PRTFile << std::setprecision(6);
    for(int32_t i=0; i<std::min(Nv, NEcho); ++i) PRTFile << std::setw(14) << v[i] << " ";
    if(Nv > NEcho) PRTFile << "... " << std::setw(14) << v[Nv-1];
    PRTFile << "\n";
}

HOST_DEVICE inline void SubTab(real *x, int32_t Nx)
{
    if(Nx >= 3){
        if(x[2] == RC(-999.9)){ // testing for equality here is dangerous
            if(x[1] == RC(-999.9)) x[1] = x[0];
            real deltax = (x[1] - x[0]) / (real)(Nx - 1);
            real x0 = x[0];
            for(int32_t i=0; i<Nx; ++i) x[i] = x0 + (real)i * deltax;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
//Copied from a different project
////////////////////////////////////////////////////////////////////////////////

#if 0
namespace math {
	//Intrinsic/optimized math functions on device, or normal ones on host.
	//http://docs.nvidia.com/cuda/cuda-math-api/group__CUDA__MATH__INTRINSIC__SINGLE.html
	//http://docs.nvidia.com/cuda/cuda-math-api/group__CUDA__MATH__SINGLE.html
	#ifdef __CUDA_ARCH__
	#define DEFINE_MATH_FUNC_1(BASE, DEV_F, HOST_F, DEV_D, HOST_D) \
		__device__ inline float BASE(const float &a) { return DEV_F(a); } \
		__device__ inline double BASE(const double &a) { return DEV_D(a); }
	#define DEFINE_MATH_FUNC_2(BASE, DEV_F, HOST_F, DEV_D, HOST_D) \
		__device__ inline float BASE(const float &a, const float &b) { return DEV_F(a, b); } \
		__device__ inline double BASE(const double &a, const double &b) { return DEV_D(a, b); }
	#else
	#define DEFINE_MATH_FUNC_1(BASE, DEV_F, HOST_F, DEV_D, HOST_D) \
		inline float BASE(const float &a) { return HOST_F(a); } \
		inline double BASE(const double &a) { return HOST_D(a); }
    #define DEFINE_MATH_FUNC_2(BASE, DEV_F, HOST_F, DEV_D, HOST_D) \
        inline float BASE(const float &a, const float &b) { return HOST_F(a, b); } \
        inline double BASE(const double &a, const double &b) { return HOST_D(a, b); }
    #endif

    #ifdef WIN32
    //These are provided as a GCC extension, but not by MSVC. But, __exp10f is a
    //device intrinsic, so an equivalent host function must exist.
    inline float internal_exp10f(float f) { return ::pow(10.0f, f); }
    inline double internal_exp10d(double d) { return ::pow(10.0, d); }
    #else
    #define internal_exp10f ::exp10
    #define internal_exp10d ::exp10
    #endif

    DEFINE_MATH_FUNC_1(abs, fabsf, std::abs, fabs, std::abs)
    DEFINE_MATH_FUNC_2(max, fmaxf, std::max, fmax, std::max)
    DEFINE_MATH_FUNC_2(min, fminf, std::min, fmin, std::min)
	DEFINE_MATH_FUNC_1(floor, floorf, std::floor, ::floor, std::floor)
	DEFINE_MATH_FUNC_1(ceil, ceilf, std::ceil, ::ceil, std::ceil)
	DEFINE_MATH_FUNC_2(fmod, fmodf, std::fmod, fmod, std::fmod)
	DEFINE_MATH_FUNC_1(sqrt, __fsqrt_rn, std::sqrt, __dsqrt_rn, std::sqrt)
	DEFINE_MATH_FUNC_1(rsqrt, __frsqrt_rn, 1.0f / std::sqrt, 1.0 / __dsqrt_rn, 1.0 / std::sqrt)
	DEFINE_MATH_FUNC_1(log2, __log2f, std::log2, log2, std::log2)
	DEFINE_MATH_FUNC_1(exp2, exp2f, std::exp2, exp2, std::exp2)
	DEFINE_MATH_FUNC_1(log10, __log10f, std::log10, log10, std::log10)
	DEFINE_MATH_FUNC_1(exp10, __exp10f, internal_exp10f, ::exp10, internal_exp10d)
	DEFINE_MATH_FUNC_2(pow, __powf, std::pow, ::pow, std::pow)
	DEFINE_MATH_FUNC_1(sin, __sinf, std::sin, ::sin, std::sin)
	DEFINE_MATH_FUNC_1(cos, __cosf, std::cos, ::cos, std::cos)
	DEFINE_MATH_FUNC_2(atan2, atan2f, std::atan2, ::atan2, std::atan2)

    //No math functions for complex on device.
    template<typename f> HOST_DEVICE inline std::complex<f> abs(std::complex<f> v){
        #ifdef __CUDA_ARCH__
        return sqrt(v.real() * v.real() + v.imag() * v.imag());
        #else
        return std::abs(v);
        #endif
    }
    template<typename f> HOST_DEVICE inline std::complex<f> sqrt(std::complex<f> v){
        #ifdef __CUDA_ARCH__
        f mag = abs(v);
        
        #else
        return std::sqrt(v);
        #endif
    }
}
#endif
