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

#ifndef _BHC_INCLUDING_COMPONENTS_
#error "Must be included from common_run.hpp!"
#endif

namespace bhc {

#if __cplusplus < 202002L
// Pre-C++20 version of bit_cast
// from https://en.cppreference.com/w/cpp/numeric/bit_cast

template<class To, class From> typename std::enable_if_t<
    sizeof(To) == sizeof(From) && std::is_trivially_copyable<From>::value
        && std::is_trivially_copyable<To>::value,
    To>
// constexpr support needs compiler magic
bit_cast(const From &src) noexcept
{
    static_assert(
        std::is_trivially_constructible<To>::value,
        "This implementation additionally requires destination type to be trivially "
        "constructible");

    To dst;
    std::memcpy(&dst, &src, sizeof(To));
    return dst;
}

#else
#include <bit>
using std::bit_cast;
#endif

inline int32_t IntFloatAdd(float v, int32_t i)
{
    return bit_cast<int32_t>(v + bit_cast<float>(i));
}

inline int64_t Int64DoubleAdd(double v, int64_t i)
{
    return bit_cast<int64_t>(v + bit_cast<double>(i));
}

//#define TECHNICALLY_UNSAFE_NONATOMIC_LOAD 1

HOST_DEVICE inline void AtomicAddReal(float *ptr, float v)
{
#ifdef __CUDA_ARCH__
    // Floating-point atomic add is natively supported on the GPU.
    atomicAdd(ptr, v);
#else
    // Have to do compare-and-swap.
    int32_t *intptr = (int32_t *)ptr;
    int32_t curint;
#ifdef __GNUC__
#ifdef TECHNICALLY_UNSAFE_NONATOMIC_LOAD
    curint = *intptr;
#else
    __atomic_load(intptr, &curint, __ATOMIC_RELAXED);
#endif
    while(!__atomic_compare_exchange_n(
        intptr, &curint, IntFloatAdd(v, curint), true, __ATOMIC_RELAXED,
        __ATOMIC_RELAXED))
        ;
#elif defined(_MSC_VER)
    int32_t prevint;
#ifdef TECHNICALLY_UNSAFE_NONATOMIC_LOAD
    curint = *intptr;
#else
    curint = InterlockedOr((LONG *)intptr, 0); // MSVC does not have a pure atomic load.
#endif
    do {
        prevint = curint;
        curint  = InterlockedCompareExchange(
            (LONG *)intptr, (LONG)IntFloatAdd(v, curint), (LONG)curint);
    } while(curint != prevint);
#else
#error "Unrecognized compiler for atomic intrinsics!"
#endif
#endif
}

HOST_DEVICE inline void AtomicAddReal(double *ptr, double v)
{
#ifdef __CUDA_ARCH__
    // Double-precision atomic add is natively supported on the GPU (compute >= 6.0).
    atomicAdd(ptr, v);
#else
    // Have to do compare-and-swap.
    int64_t *intptr = (int64_t *)ptr;
    int64_t curint;
#ifdef __GNUC__
#ifdef TECHNICALLY_UNSAFE_NONATOMIC_LOAD
    curint = *intptr;
#else
    __atomic_load(intptr, &curint, __ATOMIC_RELAXED);
#endif
    while(!__atomic_compare_exchange_n(
        intptr, &curint, Int64DoubleAdd(v, curint), true, __ATOMIC_RELAXED,
        __ATOMIC_RELAXED))
        ;
#elif defined(_MSC_VER)
    int64_t prevint;
#ifdef TECHNICALLY_UNSAFE_NONATOMIC_LOAD
    curint = *intptr;
#else
    curint = InterlockedOr64(intptr, 0);       // MSVC does not have a pure atomic load.
#endif
    do {
        prevint = curint;
        curint  = InterlockedCompareExchange64(intptr, Int64DoubleAdd(v, curint), curint);
    } while(curint != prevint);
#else
#error "Unrecognized compiler for atomic intrinsics!"
#endif
#endif
}

template<typename REAL> HOST_DEVICE inline void AtomicAddCpx(
    STD::complex<REAL> *ptr, cpx v)
{
    // This is the right way according to https://stackoverflow.com/questions/24229808/
    // getting-pointers-to-the-real-and-imaginary-parts-of-a-complex-vector-in-c
    // (though this should amount to the same thing as (real*)ptr and &((real*)ptr)[1] ).
    AtomicAddReal(&reinterpret_cast<REAL(&)[2]>(*ptr)[0], (REAL)v.real());
    AtomicAddReal(&reinterpret_cast<REAL(&)[2]>(*ptr)[1], (REAL)v.imag());
}

template<typename INT> HOST_DEVICE inline INT AtomicFetchAdd(INT *ptr, INT val)
{
#ifdef __CUDA_ARCH__
    return atomicAdd(ptr, val);
#elif defined(__GNUC__)
    return __atomic_fetch_add(ptr, val, __ATOMIC_RELAXED);
#elif defined(_MSC_VER)
    return InterlockedExchangeAdd((LONG *)ptr, (LONG)val);
#else
#error "Unrecognized compiler for atomic intrinsics!"
#endif
}

} // namespace bhc
