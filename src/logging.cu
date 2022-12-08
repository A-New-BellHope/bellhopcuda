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
#include "common.hpp"

namespace bhc {

__managed__ char gpu_log_buf[gpu_log_buf_size];
__managed__ uint32_t gpu_log_buf_pos;
uint32_t gpu_log_buf_pos_cpu;

void CudaInitLog()
{
    gpu_log_buf_pos     = 0;
    gpu_log_buf_pos_cpu = 0;
}

void CudaPostKernelLog()
{
    uint32_t gpos  = gpu_log_buf_pos & (gpu_log_buf_size - 1);
    uint32_t &cpos = gpu_log_buf_pos_cpu; // rename to shorter name
    uint32_t len   = (gpos < cpos) ? (gpos + gpu_log_buf_size - cpos) : (gpos - cpos);
    if(len == 0) return;
    char *buf = (char *)malloc(len + 1);
    buf[len]  = '\0';
    memcpy(&buf[0], &gpu_log_buf[cpos], std::min(len, gpu_log_buf_size - cpos));
    if(cpos + len > gpu_log_buf_size) {
        memcpy(&buf[gpu_log_buf_size - cpos], &gpu_log_buf[0], gpos);
    }
    GlobalLogImpl(buf);
    free(buf);
    cpos = gpos;
}

} // namespace bhc
