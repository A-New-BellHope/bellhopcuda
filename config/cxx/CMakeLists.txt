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

cmake_minimum_required(VERSION 3.9)
project(bellhopcxx LANGUAGES CXX)

include(../SetupCommon.cmake)

set(ALL_SOURCES ${COMMON_INCLUDES} ${COMMON_SOURCE} 
    ${CMAKE_SOURCE_DIR}/src/cmdline.cpp
)

if(BHC_3D_SEPARATE)
    add_executable(bellhopcxx2d ${ALL_SOURCES})
    target_compile_definitions(bellhopcxx2d PUBLIC BHC_CMDLINE=1 BHC_DIMMODE=2)
    bellhop_setup_target(bellhopcxx2d)
    add_executable(bellhopcxx3d ${ALL_SOURCES})
    target_compile_definitions(bellhopcxx3d PUBLIC BHC_CMDLINE=1 BHC_DIMMODE=3)
    bellhop_setup_target(bellhopcxx3d)
    add_executable(bellhopcxxnx2d ${ALL_SOURCES})
    target_compile_definitions(bellhopcxxnx2d PUBLIC BHC_CMDLINE=1 BHC_DIMMODE=4)
    bellhop_setup_target(bellhopcxxnx2d)
else()
    add_executable(bellhopcxx ${ALL_SOURCES})
    target_compile_definitions(bellhopcxx PUBLIC BHC_CMDLINE=1 BHC_DIMMODE=0)
    bellhop_setup_target(bellhopcxx)
endif()
