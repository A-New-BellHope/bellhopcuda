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
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_CXX_FLAGS}")

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

set(common_includes
    bhc.hpp
    math.hpp
    platform.hpp
    structs.hpp
)

set(common_source
    mode/arr.cpp
    mode/arr.hpp
    mode/eigen.cpp
    mode/eigen.hpp
    mode/field.cpp
    mode/field.hpp
    mode/fieldimpl.hpp
    mode/modemodule.hpp
    mode/ray.cpp
    mode/ray.hpp
    mode/tl.cpp
    mode/tl.hpp
    module/atten.cpp
    module/atten.hpp
    module/beaminfo.hpp
    module/botopt.hpp
    module/boundarycond.hpp
    module/boundary.hpp
    module/freq0.hpp
    module/freqvec.hpp
    module/nmedia.hpp
    module/paramsmodule.hpp
    module/rayangles.hpp
    module/rcvrbearings.hpp
    module/rcvrranges.hpp
    module/reflcoef.hpp
    module/runtype.hpp
    module/sbp.hpp
    module/ssp.hpp
    module/sxsy.hpp
    module/szrz.hpp
    module/title.hpp
    module/topopt.hpp
    util/atomics.hpp
    util/directio.hpp
    util/errors.cpp
    util/errors.hpp
    util/ldio.hpp
    util/prtfileemu.hpp
    util/timing.cpp
    util/timing.hpp
    util/unformattedio.hpp
    api.cpp
    arrivals.hpp
    boundary.hpp
    common.hpp
    common_run.hpp
    common_setup.hpp
    curves.hpp
    eigenrays.hpp
    influence.hpp
    reflect.hpp
    runtype.hpp
    ssp.hpp
    step.hpp
    trace.hpp
)

prependlist(common_includes "${CMAKE_SOURCE_DIR}/include/bhc/" ${common_includes})
prependlist(common_source "${CMAKE_SOURCE_DIR}/src/" ${common_source})

if(NOT BHC_DIM_ENABLE_2D AND NOT BHC_DIM_ENABLE_3D AND NOT BHC_DIM_ENABLE_NX2D)
    message(FATAL_ERROR "2D, 3D, and Nx2D dim modes all disabled, nothing to build!")
endif()

find_package(Threads)

function(bhc_setup_target target_name defs use_addl)
    if(BHC_USE_FLOATS)
        target_compile_definitions(${target_name} PUBLIC BHC_USE_FLOATS=1)
    endif()
    if(BHC_DEBUG)
        target_compile_definitions(${target_name} PUBLIC BHC_DEBUG=1)
    endif()
    # if(BHC_PROF AND CMAKE_COMPILER_IS_GNUCXX)
    #     target_compile_options(${target_name} PUBLIC -pg)
    #     target_link_options(${target_name} PUBLIC -pg)
    # endif()
    target_compile_definitions(${target_name} PRIVATE BHC_DLL_EXPORT=1)
    target_compile_definitions(${target_name} PUBLIC "${defs}")
    if(use_addl)
        target_compile_definitions(${target_name} PRIVATE "${addl_defs}")
        target_include_directories(${target_name} PRIVATE "${addl_includes}")
    endif()
    target_include_directories(${target_name} PUBLIC "${CMAKE_SOURCE_DIR}/include")
    target_include_directories(${target_name} PUBLIC "${CMAKE_SOURCE_DIR}/glm")
    target_link_libraries(${target_name} PUBLIC Threads::Threads)
    if(WIN32)
        set_property(TARGET ${target_name} PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded")
    endif()
endfunction()

function(bhc_create_executable target_name defs)
    add_executable(${target_name}
        $<TARGET_OBJECTS:${objlibname}>
        ${CMAKE_SOURCE_DIR}/src/cmdline.cpp
    )
    bhc_setup_target(${target_name} "${defs};BHC_CMDLINE=1" 1)
endfunction()

include(${CMAKE_SOURCE_DIR}/config/GenTemplates.cmake)

function(bhc_add_libs_exes type_name gen_extension addl_sources addl_includes addl_defs)
    set(exename "bellhop${type_name}")
    set(objlibname "${exename}objlib")
    gen_templates(${gen_extension} gen_sources)
    add_library(${objlibname} OBJECT
        ${common_includes}
        ${common_source}
        ${gen_sources}
        ${addl_sources}
    )
    set_property(TARGET ${objlibname} PROPERTY POSITION_INDEPENDENT_CODE 1)
    set(enab2d 0)
    set(enab3d 0)
    set(enabnx2d 0)
    if(BHC_DIM_ENABLE_2D)
        set(enab2d 1)
    endif()
    if(BHC_DIM_ENABLE_3D)
        set(enab3d 1)
    endif()
    if(BHC_DIM_ENABLE_NX2D)
        set(enabnx2d 1)
    endif()
    set(dim_enables "BHC_ENABLE_2D=${enab2d};BHC_ENABLE_3D=${enab3d};BHC_ENABLE_NX2D=${enabnx2d}")
    bhc_setup_target(${objlibname} "${dim_enables}" 1)
    if(BHC_LIMIT_FEATURES)
        target_compile_definitions(${objlibname} PRIVATE BHC_LIMIT_FEATURES=1)
    endif()
    add_gen_template_defs(${objlibname})
    # Targets using object library
    add_library(${exename}lib SHARED $<TARGET_OBJECTS:${objlibname}>)
    bhc_setup_target(${exename}lib "${dim_enables}" 0)
    add_library(${exename}static STATIC $<TARGET_OBJECTS:${objlibname}>)
    bhc_setup_target(${exename}static "${dim_enables}" 0)
    bhc_create_executable(${exename} "${dim_enables};BHC_DIM_ONLY=0")
    if(BHC_DIM_ENABLE_2D)
        bhc_create_executable(${exename}2d   "BHC_ENABLE_2D=1;BHC_DIM_ONLY=2")
    endif()
    if(BHC_DIM_ENABLE_3D)
        bhc_create_executable(${exename}3d   "BHC_ENABLE_3D=1;BHC_DIM_ONLY=3")
    endif()
    if(BHC_DIM_ENABLE_NX2D)
        bhc_create_executable(${exename}nx2d "BHC_ENABLE_NX2D=1;BHC_DIM_ONLY=4")
    endif()
endfunction()
