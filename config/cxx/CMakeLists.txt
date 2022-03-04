cmake_minimum_required(VERSION 3.9)
project(bellhopcxx LANGUAGES CXX)

include(../SetupCommon.cmake)

add_executable(bellhopcxx ${COMMON_INCLUDES} ${COMMON_SOURCE} 
    ${CMAKE_SOURCE_DIR}/src/cmdline.cpp
)

target_compile_definitions(bellhopcxx PUBLIC BHC_CMDLINE=1)
bellhop_setup_target(bellhopcxx)
