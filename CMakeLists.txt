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

option(BHC_DEBUG "Enable debugging features" OFF)
option(BHC_PROF  "Enable profiling features" OFF)
if(BHC_DEBUG OR BHC_PROF)
    set(BUILD_TYPE_INTERNAL RelWithDebInfo) # don't use Debug
else()
    set(BUILD_TYPE_INTERNAL Release)
endif()
set(CMAKE_BUILD_TYPE ${BUILD_TYPE_INTERNAL})
set(CMAKE_CONFIGURATION_TYPES ${CMAKE_BUILD_TYPE} CACHE STRING "" FORCE) # Only allow Visual Studio to generate the correct build type (not default to Debug)

project(bellhoptoplevel LANGUAGES NONE)

option(BHC_ENABLE_CUDA "Build CUDA version in addition to C++ version" ON)
option(BHC_ENABLE_LIB  "Build shared library version in addition to executable versions" ON)
option(BHC_3D_SEPARATE "Build separate executables for 2D/3D/Nx2D instead of selecting dimensionality with command line" OFF)
option(BHC_LIMIT_FEATURES "Limit bellhopcxx/bellhopcuda to only features supported by BELLHOP/BELLHOP3D" OFF)
option(BHC_USE_FLOATS  "Perform all floating-point arithmetic as 32-bit" OFF)

option(BHC_DIM_ENABLE_2D   "Enable 2D runs" ON)
option(BHC_DIM_ENABLE_3D   "Enable 3D runs" ON)
option(BHC_DIM_ENABLE_NX2D "Enable Nx2D runs" ON)

add_subdirectory(config)
