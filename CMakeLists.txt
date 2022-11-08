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

set(BUILD_TYPE_INTERNAL RelWithDebInfo) # Debug, RelWithDebInfo, Release

set(CMAKE_BUILD_TYPE ${BUILD_TYPE_INTERNAL})
set(CMAKE_CONFIGURATION_TYPES ${CMAKE_BUILD_TYPE} CACHE STRING "" FORCE) # Only allow Visual Studio to generate the correct build type (not default to Debug)

project(bellhoptoplevel LANGUAGES NONE)

option(BHC_ENABLE_CUDA "Build CUDA version in addition to C++ version" ON)
option(BHC_ENABLE_LIB "Build shared library version in addition to executable versions" ON)
option(BHC_3D_SEPARATE "Build separate executables for 2D/3D/Nx2D instead of selecting dimensionality with command line" OFF)
option(USE_FLOAT "Perform all floating-point arithmetic as 32-bit" OFF)

option(BHC_DIM_ENABLE_2D   "Enable 2D runs" ON)
option(BHC_DIM_ENABLE_3D   "Enable 3D runs" ON)
option(BHC_DIM_ENABLE_NX2D "Enable Nx2D runs" ON)

option(BHC_RUN_ENABLE_TL        "Enable TL runs        (Beam->RunType[0] == 'C' or 'S' or 'I')" ON)
option(BHC_RUN_ENABLE_EIGENRAYS "Enable eigenrays runs (Beam->RunType[0] == 'E')"               ON)
option(BHC_RUN_ENABLE_ARRIVALS  "Enable arrivals runs  (Beam->RunType[0] == 'A' or 'a')"        ON)

option(BHC_INFL_ENABLE_CERVENY_RAYCEN "Enable Cerveny  ray-centered influence (Beam->Type[0] == 'R')" ON)
option(BHC_INFL_ENABLE_CERVENY_CART   "Enable Cerveny  Cartesian    influence (Beam->Type[0] == 'C')" ON)
option(BHC_INFL_ENABLE_HAT_RAYCEN     "Enable hat      ray-centered influence (Beam->Type[0] == 'g')" ON)
option(BHC_INFL_ENABLE_HAT_CART       "Enable hat      Cartesian    influence (Beam->Type[0] == 'G')" ON)
option(BHC_INFL_ENABLE_GAUSS_RAYCEN   "Enable Gaussian ray-centered influence (Beam->Type[0] == 'b')" ON)
option(BHC_INFL_ENABLE_GAUSS_CART     "Enable Gaussian Cartesian    influence (Beam->Type[0] == 'B')" ON)
option(BHC_INFL_ENABLE_SGB            "Enable simple Gaussian beams influence (Beam->Type[0] == 'S')" ON)

option(BHC_SSP_ENABLE_N2LINEAR   "Enable N2-linear     1D SSP (ssp->Type == 'N')" ON)
option(BHC_SSP_ENABLE_CLINEAR    "Enable C-linear      1D SSP (ssp->Type == 'C')" ON)
option(BHC_SSP_ENABLE_CUBIC      "Enable cubic spline  1D SSP (ssp->Type == 'S')" ON)
option(BHC_SSP_ENABLE_PCHIP      "Enable PCHIP         1D SSP (ssp->Type == 'P')" ON)
option(BHC_SSP_ENABLE_QUAD       "Enable quadrilateral 2D SSP (ssp->Type == 'Q')" ON)
option(BHC_SSP_ENABLE_HEXAHEDRAL "Enable hexahedral    3D SSP (ssp->Type == 'H')" ON)
option(BHC_SSP_ENABLE_ANALYTIC   "Enable analytic      1D SSP (ssp->Type == 'A')" ON)

add_subdirectory(config)
