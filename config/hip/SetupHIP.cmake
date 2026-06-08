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

option(HIP_PRINT_REGISTERS "Print kernel register use" OFF)

set(CMAKE_HIP_STANDARD 17) # C++17
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Set default architecture if not specified
if(NOT DEFINED CMAKE_HIP_ARCHITECTURES OR CMAKE_HIP_ARCHITECTURES STREQUAL "")
    set(CMAKE_HIP_ARCHITECTURES "gfx90a")
endif()

if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set(HIP_DEBUG_FLAGS "-g")
elseif(${CMAKE_BUILD_TYPE} STREQUAL "RelWithDebInfo")
    set(HIP_DEBUG_FLAGS "-gline-tables-only")
endif()

set(HIP_EXTRA_FLAGS "-ffast-math")
if(USE_FLOATS)
    set(HIP_EXTRA_FLAGS "${HIP_EXTRA_FLAGS}")
endif()
string(STRIP "${HIP_EXTRA_FLAGS}" HIP_EXTRA_FLAGS)
message(STATUS "HIP extra flags: " ${HIP_EXTRA_FLAGS})

string(REPLACE " " ";" EXTRA_CXX_FLAGS_LIST ${EXTRA_CXX_FLAGS})
set(CMAKE_HIP_FLAGS "${HIP_DEBUG_FLAGS} ${EXTRA_CXX_FLAGS} ${HIP_EXTRA_FLAGS}")
message(STATUS "Full HIP flags: ${CMAKE_HIP_FLAGS}")
