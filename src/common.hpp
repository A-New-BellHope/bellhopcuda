#pragma once
#include <cstdio>

#ifdef __CUDACC__
#define HOST_DEVICE __host__ __device__
#include <cuda/std/complex>
#include <cuda/std/cfloat>
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

#ifdef USE_FLOATS
using real = float;
#define REAL_MAX FLT_MAX
#define REAL_EPSILON FLT_EPSILON
#define RC(a) (a##f)
//Real constant: transforms 2.0 into 2.0f. This is needed on CUDA or else
//unwanted double-precision instructions may be emitted.
#else
using real = double;
#define REAL_MAX DBL_MAX
#define REAL_EPSILON DBL_EPSILON
#define RC(a) a
#endif

using cpx = STD::complex<real>;

using vec2 = glm::vec<2, real, glm::defaultp>;
using vec3 = glm::vec<3, real, glm::defaultp>;
using vec4 = glm::vec<4, real, glm::defaultp>;
using mat2 = glm::mat<2, 2, real, glm::defaultp>;
using mat3 = glm::mat<3, 3, real, glm::defaultp>;
using mat4 = glm::mat<4, 4, real, glm::defaultp>;

#define SQ(a) ((a) * (a)) //Square


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
