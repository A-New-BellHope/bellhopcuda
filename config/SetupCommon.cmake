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

set(CMAKE_CXX_STANDARD 17) # C++17
set(CMAKE_CXX_STANDARD_REQUIRED ON) # ...is required
set(CMAKE_CXX_EXTENSTIONS OFF) # ...without compiler extensions like gnu++11
set(CMAKE_POSITION_INDEPENDENT_CODE ON) # Necessary to build shared libraries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# Set default compile flags for each platform
if(CMAKE_COMPILER_IS_GNUCXX)
    message(STATUS "GCC detected, adding compile flags")
    set(EXTRA_CXX_FLAGS "-Wall -Wextra -Wno-class-memaccess")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(STATUS "Clang detected")
    set(EXTRA_CXX_FLAGS "-Wall -Wextra")
else()
    message(STATUS "Not GCC or clang, assuming Windows format compile flags")
    # C2422: implicit conversion of double to float; BELLHOP sometimes carelessly
    # mixes doubles and floats, and we need to use exactly the same ones to match
    set(EXTRA_CXX_FLAGS "/W4 /wd4244")
endif(CMAKE_COMPILER_IS_GNUCXX)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_CXX_FLAGS}")

find_package(Threads)

include(${CMAKE_SOURCE_DIR}/config/GenTemplates.cmake)

function(bellhop_setup_target target_name)
    if(BHC_USE_FLOATS)
        target_compile_definitions(${target_name} PUBLIC BHC_USE_FLOATS=1)
    endif()
    if(BHC_LIMIT_FEATURES)
        target_compile_definitions(${target_name} PRIVATE BHC_LIMIT_FEATURES=1)
    endif()
    if(BHC_DEBUG)
        target_compile_definitions(${target_name} PRIVATE BHC_DEBUG=1)
    endif()
    add_gen_template_defs(${target_name})
    target_include_directories(${target_name} PUBLIC "${CMAKE_SOURCE_DIR}/include")
    target_include_directories(${target_name} PUBLIC "${CMAKE_SOURCE_DIR}/glm")
    target_link_libraries(${target_name} Threads::Threads)
endfunction()

function(bellhop_create_executable target_name iscuda dimmode sources defs incs)
    if(iscuda)
        set(gen_extension "cu")
    else()
        set(gen_extension "cpp")
    endif()
    gen_templates(${gen_extension} ${dimmode} gen_sources)
    add_executable(${target_name} ${sources} ${gen_sources})
    target_compile_definitions(${target_name} PUBLIC
        BHC_CMDLINE=1 BHC_DIMMODE=${dimmode} ${defs}
    )
    target_include_directories(${target_name} PRIVATE ${incs})
    bellhop_setup_target(${target_name})
endfunction()

function(bellhop_create_executables iscuda sources defs incs)
    if(iscuda)
        set(type_name "cuda")
    else()
        set(type_name "cxx")
    endif()
    if(BHC_3D_SEPARATE)
        bellhop_create_executable(bellhop${type_name}2d   ${iscuda} 2 "${sources}" "${defs}" "${incs}")
        bellhop_create_executable(bellhop${type_name}3d   ${iscuda} 3 "${sources}" "${defs}" "${incs}")
        bellhop_create_executable(bellhop${type_name}nx2d ${iscuda} 4 "${sources}" "${defs}" "${incs}")
    else()
        bellhop_create_executable(bellhop${type_name}     ${iscuda} 0 "${sources}" "${defs}" "${incs}")
    endif()
endfunction()

function(prepend OUT_VAR PREFIX) #Arguments 3, 4, etc. are items to prepend to
    set(TEMP "")
    foreach(ITEM ${ARGN})
        set(TEMP "${TEMP} ${PREFIX}${ITEM}")
    endforeach()
    set(${OUT_VAR} "${TEMP}" PARENT_SCOPE)
endfunction()

function(prependlist OUT_VAR PREFIX) #Arguments 3, 4, etc. are items to prepend to
    set(TEMP "")
    foreach(ITEM ${ARGN})
        set(TEMP "${TEMP};${PREFIX}${ITEM}")
    endforeach()
    set(${OUT_VAR} "${TEMP}" PARENT_SCOPE)
endfunction()

set(COMMON_INCLUDES
    bhc.hpp
    math.hpp
    platform.hpp
    structs.hpp
)

set(COMMON_SOURCE
    angles.hpp
    arrivals.cpp
    arrivals.hpp
    atomics.hpp
    attenuation.hpp
    beams.hpp
    bino.hpp
    boundary.hpp
    common.hpp
    curves.hpp
    eigenrays.cpp
    eigenrays.hpp
    influence.hpp
    jobs.hpp
    ldio.hpp
    logging.cpp
    logging.hpp
    prtfileemu.hpp
    raymode.hpp
    readenv.cpp
    readenv.hpp
    reflect.hpp
    run.cpp
    run.hpp
    setup.cpp
    sourcereceiver.hpp
    ssp.cpp
    ssp.hpp
    step.hpp
    tlmode.cpp
    tlmode.hpp
    trace.hpp
)

prependlist(COMMON_INCLUDES "${CMAKE_SOURCE_DIR}/include/bhc/" ${COMMON_INCLUDES})
prependlist(COMMON_SOURCE "${CMAKE_SOURCE_DIR}/src/" ${COMMON_SOURCE})
