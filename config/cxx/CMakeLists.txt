cmake_minimum_required(VERSION 3.9)
project(bellhopcxx LANGUAGES CXX)

include(../SetupCommon.cmake)

add_executable(bellhopcxx ${COMMON_SOURCE} ${CMAKE_SOURCE_DIR}/src/tl_cxx.cpp)

bellhop_setup_target(bellhopcxx)
