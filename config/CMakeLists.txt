add_subdirectory(cxx)

if(NOT(APPLE OR DEFINED ENV{BELLHOP_NO_CUDA}))
    add_subdirectory(cuda)
endif()
