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

#include "common.hpp"

#include <cctype>
#include <vector>
#include <string>
#include <locale>
#include <algorithm>
//#include <cfenv>
#include <exception>

// More includes below.

namespace bhc {

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

HOST_DEVICE inline void PrintMatrix(const mat2x2 &m, const char *label)
{
    printf(
        "%s: /%12.7e %12.7e\\\n       \\%12.7e %12.7e/\n", label, m[0][0], m[1][0],
        m[0][1], m[1][1]);
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

inline bool startswith(const std::string &source, const std::string &target)
{
    return source.rfind(target, 0) == 0;
}

inline bool endswith(const std::string &source, const std::string &target)
{
    size_t l = source.length() - target.length();
    return source.find(target, l) == l;
}

} // namespace bhc

#define _BHC_INCLUDING_COMPONENTS_ 1
#include "util/ldio.hpp"
#include "util/directio.hpp"
#include "util/unformattedio.hpp"
#undef _BHC_INCLUDING_COMPONENTS_

namespace bhc {

////////////////////////////////////////////////////////////////////////////////
// CUDA memory
////////////////////////////////////////////////////////////////////////////////

template<bool O3D, typename T> inline void trackdeallocate(
    const bhcParams<O3D> &params, T *&ptr)
{
    if(ptr == nullptr) return;
    // Size stored two 64-bit words before returned pointer. 16 byte aligned.
    uint64_t *ptr2 = (uint64_t *)ptr;
    ptr2 -= 2;
    GetInternal(params)->usedMemory -= *ptr2;
#ifdef BHC_BUILD_CUDA
    checkCudaErrors(cudaFree(ptr2));
#else
    free(ptr2);
#endif
    ptr = nullptr;
}

template<bool O3D, typename T> inline void trackallocate(
    const bhcParams<O3D> &params, const char *description, T *&ptr, size_t n = 1)
{
    if(ptr != nullptr) trackdeallocate(params, ptr);
    uint64_t *ptr2;
    uint64_t s  = ((n * sizeof(T)) + 15ull) & ~15ull; // Round up to 16 byte aligned
    uint64_t s2 = s + 16ull; // Total size to allocate, including size info
    if(GetInternal(params)->usedMemory + s2 > GetInternal(params)->maxMemory) {
        EXTERR(
            "Insufficient memory to allocate %s, need more than %" PRIu64 " MiB",
            description, (GetInternal(params)->usedMemory + s2) / (1024ull * 1024ull));
    }
#ifdef BHC_BUILD_CUDA
    checkCudaErrors(cudaMallocManaged(&ptr2, s2));
#else
    ptr2 = (uint64_t *)malloc(s2);
#endif
    *ptr2 = s2;
    GetInternal(params)->usedMemory += s2;
    ptr = (T *)(ptr2 + 2);
#ifdef BHC_DEBUG
    // Debugging: Fill memory with garbage data to help detect uninitialized vars
    memset(ptr, 0xFE, s);
#endif
}

////////////////////////////////////////////////////////////////////////////////
// Vector input related
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
template<typename REAL> HOST_DEVICE inline bool monotonic(
    REAL *arr, size_t n, bool repeatOK = false)
{
    CHECK_REAL_T();
    if(n < 2) return true;
    for(size_t i = 0; i < n - 1; ++i) {
        REAL next = arr[i + 1];
        REAL cur  = arr[i];
        if(next < cur || (next == cur && !repeatOK)) return false;
    }
    return true;
}
template<typename VEC2> HOST_DEVICE inline bool monotonic(
    VEC2 *arr, size_t n, const int32_t dim)
{
    if(n < 2) return true;
    for(size_t i = 0; i < n - 1; ++i) {
        if(arr[i + 1][dim] <= arr[i][dim]) return false;
    }
    return true;
}
template<typename REAL> HOST_DEVICE inline bool monotonic(
    REAL *arr, size_t n, const int32_t stridereals, const int32_t offset)
{
    CHECK_REAL_T();
    if(n < 2) return true;
    for(size_t i = 0; i < n - 1; ++i) {
        if(arr[(i + 1) * stridereals + offset] <= arr[i * stridereals + offset])
            return false;
    }
    return true;
}

/**
 * mbp: full 360-degree sweep? remove duplicate angle/beam
 */
template<typename REAL> HOST_DEVICE inline void CheckFix360Sweep(
    const REAL *angles, int32_t &n, bool &removed)
{
    CHECK_REAL_T();
    // LP: Changed from (the FORTRAN equivalent of) REAL_MINPOS, see Fortran
    // version readme.
    if(n > 1
       && STD::abs(STD::fmod(angles[n - 1] - angles[0], FL(360.0)))
           < FL(10.0) * spacing(RL(360.0))) {
        --n;
        removed = true;
    }
}

template<typename REAL> inline void ToMeters(REAL *&x, int32_t &Nx)
{
    for(int32_t i = 0; i < Nx; ++i) x[i] *= FL(1000.0);
}

/**
 * If x[2] == -999.9 then subtabulation is performed
 * i.e., a vector is generated with Nx points in [x[0], x[1]]
 * If x[1] == -999.9 then x[0] is repeated into x[1]
 */
template<typename REAL> inline void SubTab(REAL *x, int32_t Nx)
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

/**
 * Reverse the operation of SubTab. Write a vector's length and contents to
 * the list-directed output file, but if the contents are uniformly spaced,
 * only write the first and last values and the slash to stop reading.
 */
template<typename REAL> inline void UnSubTab(
    LDOFile &f, REAL *x, int32_t Nx, const char *nLabel, const char *xLabel,
    REAL scale = (REAL)(1), int32_t stridereals = 1)
{
    bool doSubTab = true;
    REAL firstx = (REAL)(0), lastx = (REAL)(0);
    if(Nx <= 2) {
        doSubTab = false;
    } else {
        firstx      = x[0];
        lastx       = x[stridereals * (Nx - 1)];
        REAL deltax = (lastx - firstx) / (REAL)(Nx - 1);
        REAL curx   = firstx;
        REAL thresh = (REAL)(3 * Nx)
            * spacing(std::max(std::abs(firstx), std::abs(lastx)));
        thresh = std::max(thresh, (REAL)(1e-7));
        for(int32_t i = 1; i < Nx - 1; ++i) {
            curx += deltax;
            if(abs(x[stridereals * i] - curx) > thresh) {
                doSubTab = false;
                break;
            }
        }
    }
    f << Nx;
    if(nLabel != nullptr) {
        f.write("! ");
        f.write(nLabel);
    }
    f << '\n';
    if(doSubTab) {
        f << (firstx * scale) << (lastx * scale) << '/' << ' ';
    } else {
        for(int32_t i = 0; i < Nx; ++i) { f << (x[stridereals * i] * scale); }
    }
    if(xLabel != nullptr) {
        f.write("! ");
        f.write(xLabel);
    }
    f << '\n';
}

/**
 * Read a vector x
 */
template<bool O3D, typename REAL> inline void ReadVector(
    bhcParams<O3D> &params, REAL *&x, int32_t &Nx, LDIFile &ENVFile,
    const char *Description)
{
    LIST(ENVFile);
    ENVFile.Read(Nx);
    trackallocate(params, Description, x, bhc::max(3, Nx));
    x[1] = FL(-999.9);
    x[2] = FL(-999.9);
    LIST(ENVFile);
    ENVFile.Read(x, Nx);
    SubTab(x, Nx);
    Sort(x, Nx);
}

template<bool O3D, typename REAL> inline void ValidateVector(
    bhcParams<O3D> &params, REAL *&x, int32_t &Nx, const char *Description,
    bool repeatOK = false)
{
    if(Nx <= 0) { EXTERR("ValidateVector: Number of %s must be positive", Description); }
    if(!monotonic(x, Nx, repeatOK)) {
        EXTERR("ValidateVector: %s are not monotonically increasing", Description);
    }
}

template<typename REAL> inline void EchoVector(
    REAL *v, int32_t Nv, PrintFileEmu &PRTFile, int32_t NEcho = 10,
    const char *ExtraSpaces = "", REAL multiplier = RL(1.0), int32_t stridereals = 1,
    int32_t offset = 0)
{
    CHECK_REAL_T();
    PRTFile << std::setprecision(6) << ExtraSpaces;
    for(int32_t i = 0, r = 0; i < bhc::min(Nv, NEcho); ++i) {
        PRTFile << std::setw(14) << (multiplier * v[i * stridereals + offset]) << " ";
        ++r;
        if(r == 5) {
            r = 0;
            PRTFile << "\n" << ExtraSpaces;
        }
    }
    if(Nv > NEcho)
        PRTFile << "... " << std::setw(14)
                << (multiplier * v[(Nv - 1) * stridereals + offset]);
    PRTFile << "\n";
    /*
    PRTFile << std::setprecision(12);
    for(int32_t i=0; i<Nv; ++i) PRTFile << std::setw(20) << v[i] << "\n";
    PRTFile << "\n";
    */
}

/**
 * LP: Echo vector with description
 * Description is something like 'receiver ranges'
 * Units       is something like 'km'
 */
template<bool O3D, typename REAL> inline void EchoVectorWDescr(
    bhcParams<O3D> &params, REAL *&x, int32_t &Nx, REAL multiplier,
    const char *Description, const char *Units)
{
    PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;

    PRTFile << "\n_______________________________________________________________________"
               "___\n\n";
    PRTFile << "   Number of " << Description << " = " << Nx << "\n";
    PRTFile << "   " << Description << " (" << Units << ")\n";
    EchoVector(x, Nx, PRTFile, 10, "   ", multiplier);
    PRTFile << "\n";
}

} // namespace bhc
