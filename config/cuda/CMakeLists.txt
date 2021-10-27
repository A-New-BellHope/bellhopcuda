project(bellhopcuda LANGUAGES CXX CUDA)

include(../SetupCommon.cmake)

include(SetupCUDA.cmake)



target_compile_definitions(TODO PUBLIC BUILD_CUDA=1)
