cmake_minimum_required(VERSION 3.9)
project(bellhoptoplevel LANGUAGES NONE)

option(USE_FLOAT "Perform all floating-point arithmetic as 32-bit" OFF)

add_subdirectory(config)
