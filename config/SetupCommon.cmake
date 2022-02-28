
set(CMAKE_CXX_STANDARD 14) # C++14
set(CMAKE_CXX_STANDARD_REQUIRED ON) # ...is required
set(CMAKE_CXX_EXTENSTIONS OFF) # ...without compiler extensions like gnu++11
set(CMAKE_POSITION_INDEPENDENT_CODE ON) # Necessary to build shared libraries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# Set default compile flags for each platform
if(CMAKE_COMPILER_IS_GNUCXX)
    message(STATUS "GCC detected, adding compile flags")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra \
    -Wno-sign-compare -Wno-unused-parameter -Wno-class-memaccess")
else()
    message(STATUS "Not GCC, assuming Windows format compile flags")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W4")
endif(CMAKE_COMPILER_IS_GNUCXX)

find_package(Threads)

function(bellhop_setup_target target_name)
    if(USE_FLOAT)
        target_compile_definitions(${target_name} PUBLIC USE_FLOATS=1)
    endif()
    target_include_directories(${target_name} PUBLIC "${CMAKE_SOURCE_DIR}/glm")
    target_link_libraries(${target_name} Threads::Threads)
endfunction()

function(prepend OUT_VAR PREFIX) #Arguments 3, 4, etc. are items to prepend to
    set(TEMP "")
    foreach(ITEM ${ARGN})
        set(TEMP "${TEMP} ${PREFIX}${ITEM}")
    endforeach()
    set(${OUT_VAR} "${TEMP}" PARENT_SCOPE)
endfunction()

function(prependlist OUT_VAR PREFIX) #Arguments 3, 4, etc. are items to prepend to
    set(TEMP "")
    foreach(ITEM ${ARGN})
        set(TEMP "${TEMP};${PREFIX}${ITEM}")
    endforeach()
    set(${OUT_VAR} "${TEMP}" PARENT_SCOPE)
endfunction()

set(COMMON_SOURCE
    angles.hpp
    arrivals.cpp
    arrivals.hpp
    atomics.hpp
    attenuation.hpp
    beams.hpp
    boundary.hpp
    common.hpp
    curves.hpp
    eigenrays.cpp
    eigenrays.hpp
    influence.hpp
    ldio.hpp
    main.hpp
    output.cpp
    output.hpp
    readenv.cpp
    readenv.hpp
    refcoef.hpp
    setup.cpp
    setup.hpp
    sourcereceiver.hpp
    ssp.cpp
    ssp.hpp
    step.hpp
    trace.hpp
)

prependlist(COMMON_SOURCE "${CMAKE_SOURCE_DIR}/src/" ${COMMON_SOURCE})
