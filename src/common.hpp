#pragma once

#include <bhc/bhc.hpp>

#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>

#ifdef BUILD_CUDA
#include "UtilsCUDA.cuh"
#endif

////////////////////////////////////////////////////////////////////////////////
//Assertions and debug
////////////////////////////////////////////////////////////////////////////////

#define NULLSTATEMENT ((void)0)
#define REQUIRESEMICOLON do{NULLSTATEMENT;} while(false)

/**
 * Returns a pointer to only the last portion of the source filename.
 * This works perfectly correctly on GPU, but it consumes many registers,
 * which are often evaluated once at the beginning and left in registers
 * through the whole kernel.
*/
/*HOST_DEVICE*/ inline const char *SOURCE_FILENAME(const char *file){
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

#ifdef __CUDA_ARCH__
#define bail __trap
#define BASSERT_STR(x) #x
#define BASSERT_XSTR(x) BASSERT_STR(x)
#define BASSERT(statement) \
if(__builtin_expect(!(statement), 0)) { \
	printf("Assertion " #statement " failed line " BASSERT_XSTR(__LINE__) "!\n"); \
	__trap(); \
} REQUIRESEMICOLON
#else
#define bail std::abort
#define BASSERT(statement) \
if(!(statement)){ \
	std::cout << "FATAL: Assertion \"" #statement "\" failed in " \
		<< SOURCE_FILENAME(__FILE__) << " line " << __LINE__ << "\n"; \
	std::abort(); \
} REQUIRESEMICOLON
#endif

////////////////////////////////////////////////////////////////////////////////
//Misc real math
////////////////////////////////////////////////////////////////////////////////

#define SQ(a) ((a) * (a)) //Square
#define CUBE(a) ((a) * (a) * (a))
constexpr real RadDeg = RL(180.0) / REAL_PI;
constexpr real DegRad = REAL_PI / RL(180.0);

template<typename REAL> HOST_DEVICE inline REAL spacing(REAL v){
	return STD::abs(STD::nextafter(v, (REAL)(0.0f)) - v);
}

namespace math {
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
}

////////////////////////////////////////////////////////////////////////////////
//String manipulation
////////////////////////////////////////////////////////////////////////////////

#include <cctype>
#include <cstring>
#include <string>
#include <locale>

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

#include <algorithm>

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
	T* ret;
	#ifdef BUILD_CUDA
		checkCudaErrors(cudaMallocManaged(&ret, n*sizeof(T)));
	#else
		ret = new T[n];
	#endif
	// Debugging: Fill memory with garbage data to help detect uninitialized vars
	memset(ret, 0xFE, n*sizeof(T));
	return ret;
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
		if(arr[(i+1)*stridereals+offset] <= arr[i*stridereals+offset]) return false;
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

template<typename REAL> inline void EchoVector(REAL *v, int32_t Nv, std::ofstream &PRTFile)
{
    constexpr int32_t NEcho = 10;
    PRTFile << std::setprecision(6);
    for(int32_t i=0, r=0; i<math::min(Nv, NEcho); ++i){
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

//#include <cfenv>

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

#include <chrono>

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

////////////////////////////////////////////////////////////////////////////////
//Other components
////////////////////////////////////////////////////////////////////////////////

#define _BHC_INCLUDING_COMPONENTS_ 1
#include "ldio.hpp"
#include "bino.hpp"
#include "atomics.hpp"
#undef _BHC_INCLUDING_COMPONENTS_
