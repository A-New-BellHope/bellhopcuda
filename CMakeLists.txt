cmake_minimum_required(VERSION 3.9)

set(BUILD_TYPE_INTERNAL RelWithDebInfo) # Debug, RelWithDebInfo, Release

set(CMAKE_BUILD_TYPE ${BUILD_TYPE_INTERNAL})
set(CMAKE_CONFIGURATION_TYPES ${CMAKE_BUILD_TYPE} CACHE STRING "" FORCE) # Only allow Visual Studio to generate the correct build type (not default to Debug)

project(bellhoptoplevel LANGUAGES NONE)

option(USE_FLOAT "Perform all floating-point arithmetic as 32-bit" OFF)

add_subdirectory(config)
