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

project(bellhopcuda LANGUAGES CXX)

include(../SetupCommon.cmake)
include(SetupHIP.cmake)

set(addl_sources
    ${CMAKE_SOURCE_DIR}/src/util/UtilsCUDA.cuh
    ${CMAKE_SOURCE_DIR}/src/util/cuda_to_hip.h
)

# Find HIP include directory
find_path(HIP_INCLUDE_DIR hip/hip_runtime.h
    HINTS
        $ENV{ROCM_PATH}/include
        $ENV{HIP_PATH}/include
        /opt/rocm/include
)

if(NOT HIP_INCLUDE_DIR)
    message(FATAL_ERROR "Could not find HIP include directory")
endif()

bhc_add_libs_exes(cuda cu "${addl_sources}" "${HIP_INCLUDE_DIR}" "BHC_BUILD_CUDA=1;USE_HIP=1")
