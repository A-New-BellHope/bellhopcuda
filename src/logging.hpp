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

#ifndef _BHC_INCLUDING_COMPONENTS_
#error "Must be included from common.hpp!"
#endif

#include <cstdarg>

namespace bhc {

#ifdef BHC_BUILD_CUDA
// Must be power of 2 so overflows properly
constexpr uint32_t gpu_log_buf_size_bits = 20;
constexpr uint32_t gpu_log_buf_size = 1 << gpu_log_buf_size_bits;
extern __managed__ char gpu_log_buf[gpu_log_buf_size]; // Circular buffer
extern __managed__ uint32_t gpu_log_buf_pos;
extern uint32_t gpu_log_buf_pos_cpu;

__device__ inline int strlen_d(const char *d){
    int ret = 0;
    while(d[ret]) ++ret;
    return ret;
}

#endif

extern void (*external_global_log)(const char *message);

inline void GlobalLogImpl(const char *message){
    if(external_global_log != nullptr){
        external_global_log(message);
    }else{
        printf("%s", message);
    }
}

/**
 * Replacement for printf which can be redirected in Unreal (including on GPU).
 * TODO does not actually print the formatted contents on GPU, just the format
 * string, due to vsnprintf etc. not being available on the GPU and the
 * complexity of using a third-party library for this.
 */
HOST_DEVICE inline void GlobalLog(const char *message, ...){
    #ifdef __CUDA_ARCH__
    uint32_t outlen = strlen_d(message);
    uint32_t pos = atomicAdd(&gpu_log_buf_pos, outlen) & (gpu_log_buf_size - 1);
    memcpy(&gpu_log_buf[gpu_log_buf_pos], &message[0],
        ::min((uint32_t)outlen, (uint32_t)(gpu_log_buf_size - pos)));
    if(pos + outlen > gpu_log_buf_size){
        memcpy(&gpu_log_buf[0], &message[gpu_log_buf_size - pos],
            pos + outlen - gpu_log_buf_size);
    }
    #else
    constexpr uint32_t maxbufsize = 1024;
    char buf[maxbufsize];
    va_list argp;
    va_start(argp, message);
    vsnprintf(buf, maxbufsize, message, argp);
    va_end(argp);
    GlobalLogImpl(buf);
    #endif
}

inline void InitLog(void (*outputCallback)(const char *message)){
    external_global_log = outputCallback;
    #ifdef BHC_BUILD_CUDA
    gpu_log_buf_pos = 0;
    gpu_log_buf_pos_cpu = 0;
    #endif
}

#ifdef BHC_BUILD_CUDA
void CudaPostKernelLog();
#endif

}
