add_subdirectory(cxx)

#if(MSVC) # this is supposed to work, but its not recognizing
if(NOT (CMAKE_COMPILER_IS_GNUCXX OR APPLE))
    message(STATUS "Adding build for Windows DLL")
    add_subdirectory(winlib)
endif()

if(NOT(APPLE OR DEFINED ENV{BELLHOP_NO_CUDA}))
    add_subdirectory(cuda)
endif()
