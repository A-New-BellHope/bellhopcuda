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
project(bellhopcxxlib LANGUAGES CXX)

include(../SetupCommon.cmake)

add_library(bellhopcxxlib SHARED ${COMMON_SOURCE}
    #no additional source files
)
#LP: If this is being used in a parent CMake project, just
#target_link_libraries(parentproject bellhopcxxlib) should work and get the
#includes. If this is being used with another build system, e.g. Unreal, just
#copy the include and glm directories into your build system and set them as
#public include paths.
#include(GenerateExportHeader)
#generate_export_header(bellhopcxxlib)

target_compile_definitions(bellhopcxxlib PUBLIC BHC_EXPORTS=1 BHC_CPU=1)
bellhop_setup_target(bellhopcxxlib)
