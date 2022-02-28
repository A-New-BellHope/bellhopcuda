cmake_minimum_required(VERSION 3.9)
project(bellhopcxxlib LANGUAGES CXX)

include(../SetupCommon.cmake)

add_library(bellhopcxxlib SHARED ${COMMON_SOURCE} ${CMAKE_SOURCE_DIR}/src/trace_lib.cpp)
include(GenerateExportHeader)
generate_export_header(bellhopcxxlib)

bellhop_setup_target(bellhopcxxlib)
