cmake_policy(SET CMP0104 OLD)
project(bellhopcuda LANGUAGES CXX CUDA)

include(../SetupCommon.cmake)
include(SetupCUDA.cmake)

add_executable(bellhopcuda ${COMMON_SOURCE} 
    ${CMAKE_SOURCE_DIR}/src/UtilsCUDA.cuh
    ${CMAKE_SOURCE_DIR}/src/tl_cuda.cu
)

target_compile_definitions(bellhopcuda PUBLIC BUILD_CUDA=1)
bellhop_setup_target(bellhopcuda)
