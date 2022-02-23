cmake_minimum_required(VERSION 3.9)
project(bellhopcxx LANGUAGES CXX)

include(../SetupCommon.cmake)

add_library(bellhopcxxdll SHARED ${COMMON_SOURCE} ${CMAKE_SOURCE_DIR}/src/trace_dll.cpp)
include(GenerateExportHeader)
generate_export_header(bellhopcxxdll)

bellhop_setup_target(bellhopcxxdll)
