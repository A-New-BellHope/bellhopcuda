# bellhopcxx / bellhopcuda - C++/CUDA port of BELLHOP underwater acoustics simulator
# Copyright (C) 2021-2022 The Regents of the University of California
# c/o Jules Jaffe team at SIO / UCSD, jjaffe@ucsd.edu
# Based on BELLHOP, which is Copyright (C) 1983-2020 Michael B. Porter
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

cmake_policy(SET CMP0104 OLD)
project(bellhopcuda LANGUAGES CXX CUDA)

include(../SetupCommon.cmake)
include(SetupCUDA.cmake)

set(ALL_SOURCES ${COMMON_INCLUDES} ${COMMON_SOURCE} 
    ${CMAKE_SOURCE_DIR}/src/cmdline.cpp
    ${CMAKE_SOURCE_DIR}/src/UtilsCUDA.cuh
    ${CMAKE_SOURCE_DIR}/src/logging.cu
)

bellhop_create_executables(1 "${ALL_SOURCES}"
    "BHC_BUILD_CUDA=1"
    "${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}"
)
