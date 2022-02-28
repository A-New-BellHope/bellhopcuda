add_subdirectory(cxx)
add_subdirectory(cxxlib)

if(NOT(APPLE OR DEFINED ENV{BELLHOP_NO_CUDA}))
    add_subdirectory(cuda)
endif()
