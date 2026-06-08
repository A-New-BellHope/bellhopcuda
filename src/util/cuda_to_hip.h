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

// CUDA-to-HIP compatibility header for bellhopcuda
// On ROCm builds, this header aliases CUDA API symbols to their HIP equivalents.
// On NVIDIA builds, this is a no-op include of the CUDA runtime.

#if defined(USE_HIP) || defined(__HIP_PLATFORM_AMD__)

#include <hip/hip_runtime.h>

// Error types and success codes
#define cudaError_t                 hipError_t
#define cudaSuccess                 hipSuccess

// Error handling
#define cudaGetErrorName            hipGetErrorName
#define cudaGetErrorString          hipGetErrorString
#define cudaGetLastError            hipGetLastError
#define cudaPeekAtLastError         hipPeekAtLastError

// Device management
#define cudaGetDeviceCount          hipGetDeviceCount
#define cudaGetDeviceProperties     hipGetDeviceProperties
#define cudaSetDevice               hipSetDevice
#define cudaDeviceSynchronize       hipDeviceSynchronize
#define cudaDeviceReset             hipDeviceReset
#define cudaDeviceProp              hipDeviceProp_t

// Memory management
#define cudaMalloc                  hipMalloc
#define cudaMallocManaged           hipMallocManaged
#define cudaFree                    hipFree
#define cudaMemcpy                  hipMemcpy
#define cudaMemcpyHostToDevice      hipMemcpyHostToDevice
#define cudaMemcpyDeviceToHost      hipMemcpyDeviceToHost

// Kernel launch bounds are already defined by HIP

#else

#include <cuda_runtime.h>

#endif
