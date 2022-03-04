cmake_minimum_required(VERSION 3.9)
project(bellhopcxxlib LANGUAGES CXX)

include(../SetupCommon.cmake)

add_library(bellhopcxxlib SHARED ${COMMON_SOURCE}
    #no additional source files
)
#LP: If this is being used in a parent CMake project, just
#target_link_libraries(parentproject bellhopcxxlib) should work and get the
#includes. If this is being used with another build system, e.g. Unreal, just
#copy the include and glm directories into your build system and set them as
#public include paths.
#include(GenerateExportHeader)
#generate_export_header(bellhopcxxlib)

target_compile_definitions(bellhopcxxlib PUBLIC BHC_EXPORTS=1)
bellhop_setup_target(bellhopcxxlib)
