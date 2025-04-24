# bellhopcxx / bellhopcuda - C++/CUDA port of BELLHOP / BELLHOP3D underwater acoustics simulator
# Copyright (C) 2021-2023 The Regents of the University of California
# Marine Physical Lab at Scripps Oceanography, c/o Jules Jaffe, jjaffe@ucsd.edu
# Based on BELLHOP / BELLHOP3D, which is Copyright (C) 1983-2022 Michael B. Porter
# 
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <https://www.gnu.org/licenses/>.

option(CUDA_PRINT_REGISTERS "Print kernel register use" OFF)
option(CUDA_DISASSEMBLY "Save temp outputs for disassembly" OFF)
##set(CUDA_ARCH_OVERRIDE "" CACHE STRING "Compile for this GPU architecture (e.g. 86)")

set(CMAKE_CUDA_STANDARD 17) # C++17
#set(CMAKE_CUDA_SEPARABLE_COMPILATION ON) # Issues on Windows
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CUDA_ARCHITECTURES native)

if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
	set(NVCC_DEBUG_FLAGS "-g -G")
elseif(${CMAKE_BUILD_TYPE} STREQUAL "RelWithDebInfo")
	set(NVCC_DEBUG_FLAGS "-lineinfo")
endif()
if(CUDA_PRINT_REGISTERS)
	set(NVCC_DEBUG_FLAGS "${NVCC_DEBUG_FLAGS} --ptxas-options=-v")
endif()
string(STRIP "${NVCC_DEBUG_FLAGS}" NVCC_DEBUG_FLAGS)
if(CUDA_DISASSEMBLY)
	set(CUDA_EXTRA_FLAGS "--keep --source-in-ptx")
endif()

set(CUDAFE_WARNINGS
    local_variable_hidden
    #variable_hides_entity
    decl_hides_catch_parameter
    decl_hides_function_parameter
    decl_hides_template_parameter
    for_init_hides_declaration
    declaration_hides_for_init
)
prepend(CUDAFE_WARNINGS_POST "-Xcudafe --diag_warning=" ${CUDAFE_WARNINGS})
set(CUDAFE_FLAGS "${CUDAFE_WARNINGS_POST} -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored")
string(STRIP "${CUDAFE_FLAGS}" CUDAFE_FLAGS)

set(CUDA_EXTRA_FLAGS "${CUDA_EXTRA_FLAGS} ${CUDAFE_FLAGS} -ftz=true --expt-relaxed-constexpr")
if(USE_FLOATS)
    set(CUDA_EXTRA_FLAGS "${CUDA_EXTRA_FLAGS} -use_fast_math -prec-div=false -prec-sqrt=false")
endif()
string(STRIP "${CUDA_EXTRA_FLAGS}" CUDA_EXTRA_FLAGS)
message(STATUS "CUDA extra flags: " ${CUDA_EXTRA_FLAGS})

string(REPLACE " " ";" EXTRA_CXX_FLAGS_LIST ${EXTRA_CXX_FLAGS})
prepend(CXX_TO_CUDA_FLAGS "-Xcompiler=" ${EXTRA_CXX_FLAGS_LIST})
set(CMAKE_CUDA_FLAGS "${NVCC_DEBUG_FLAGS} ${CXX_TO_CUDA_FLAGS} ${CUDA_EXTRA_FLAGS}")
message(STATUS "Full CUDA flags: ${CMAKE_CUDA_FLAGS}")
