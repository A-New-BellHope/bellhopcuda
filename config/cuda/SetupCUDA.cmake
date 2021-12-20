option(CUDA_PRINT_REGISTERS "Print kernel register use" OFF)
option(CUDA_DISASSEMBLY "Save temp outputs for disassembly" OFF)
set(CUDA_ARCH_OVERRIDE "" CACHE STRING "Compile for this GPU architecture (e.g. 86)")

set(CMAKE_CUDA_STANDARD 14) # C++14

if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
	set(NVCC_DEBUG_FLAGS "-g -G")
elseif(${CMAKE_BUILD_TYPE} STREQUAL "RelWithDebInfo")
	set(NVCC_DEBUG_FLAGS "-lineinfo")
endif()
if(CUDA_PRINT_REGISTERS)
	set(NVCC_DEBUG_FLAGS "${NVCC_DEBUG_FLAGS} --ptxas-options=-v")
endif()
string(STRIP "${NVCC_DEBUG_FLAGS}" NVCC_DEBUG_FLAGS)
if(CUDA_DISASSEMBLY)
	set(CUDA_EXTRA_FLAGS "--keep --source-in-ptx")
endif()

set(CUDAFE_WARNINGS
    local_variable_hidden
    #variable_hides_entity
    decl_hides_catch_parameter
    decl_hides_function_parameter
    decl_hides_template_parameter
    for_init_hides_declaration
    declaration_hides_for_init
)
prepend(CUDAFE_WARNINGS_POST "-Xcudafe --diag_warning=" ${CUDAFE_WARNINGS})
set(CUDAFE_FLAGS "${CUDAFE_WARNINGS_POST} -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored")
string(STRIP "${CUDAFE_FLAGS}" CUDAFE_FLAGS)

set(CUDA_EXTRA_FLAGS "${CUDA_EXTRA_FLAGS} ${CUDAFE_FLAGS} -ftz=true --expt-relaxed-constexpr")
if(USE_FLOATS)
    set(CUDA_EXTRA_FLAGS "${CUDA_EXTRA_FLAGS} -use_fast_math -prec-div=false -prec-sqrt=false")
endif()
string(STRIP "${CUDA_EXTRA_FLAGS}" CUDA_EXTRA_FLAGS)
message(STATUS "CUDA extra flags: " ${CUDA_EXTRA_FLAGS})

# Determine the GPU compute capabilities / gencodes to compile for
function(get_compute_for_gpu OUT_VAR GPU_NAME)
	set(GPU_DATABASE "30:GTX 770,GTX 760,GT 740,GTX 690,GTX 680,GTX 670,GTX 660 Ti,GTX 660,GTX 650 Ti BOOST,GTX 650 Ti,GTX 650;35:GTX Titan Z,GTX Titan Black,GTX Titan,GTX 780 Ti,GTX 780;50:GTX 750 Ti,GTX 750;52:GTX Titan X,GTX 980 Ti,GTX 980,GTX 970,GTX 960,GTX 950;61:TITAN Xp,Titan X,GTX 1080 Ti,GTX 1080,GTX 1070 Ti,GTX 1070,GTX 1070 with Max-Q Design,GTX 1060,GTX 1050 Ti,GTX 1050;70:TITAN V;75:TITAN RTX,RTX 2080 Ti,RTX 2080 Super,RTX 2080,RTX 2070 Super,RTX 2070,RTX 2060 Super,RTX 2060,GTX 1660 Ti,GTX 1660 Super,GTX 1660,GTX 1650 Super,GTX 1650;86:RTX 3090,RTX 3080 Ti,RTX 3080,RTX 3070 Ti,RTX 3070,RTX 3060 Ti,RTX 3060")
    foreach(COMPUTE_ROW ${GPU_DATABASE})
        string(REGEX REPLACE ":.*" "" COMPUTE_CAPABILITY ${COMPUTE_ROW})
        string(REGEX REPLACE ".*:" "" GPU_LIST ${COMPUTE_ROW})
        string(REGEX REPLACE "," ";" GPU_LIST ${GPU_LIST})
        if(${GPU_NAME} IN_LIST GPU_LIST)
            set(${OUT_VAR} ${COMPUTE_CAPABILITY} PARENT_SCOPE)
            return()
        endif()
    endforeach()
    set(${OUT_VAR} "Unknown" PARENT_SCOPE)
    message(SEND_ERROR "CUDA compute capability of GPU " ${GPU_NAME} " unknown!")
endfunction()

function(get_gencode_args OUT_VAR)
    set(TEMP "")
    if(CUDA_ARCH_OVERRIDE)
        string(REGEX REPLACE "," ";" CUDA_ARCH_OVERRIDE ${CUDA_ARCH_OVERRIDE})
        foreach(ARCH ${CUDA_ARCH_OVERRIDE})
            if(ARCH MATCHES "^[0-9]+$" AND ARCH GREATER_EQUAL 30 AND ARCH LESS_EQUAL 75)
                message(STATUS "GPU architecture override: ${ARCH}")
                set(TEMP "${TEMP} -gencode arch=compute_${ARCH},code=sm_${ARCH}")
            else()
                message(FATAL_ERROR "Invalid GPU architecure override: ${ARCH}!")
            endif()
        endforeach()
    else()
    	if(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    		set(NVIDIA_SMI_COMMAND "nvidia-smi.exe")
            set(NVIDIA_SMI_FAILURE_MSG "Could not find nvidia-smi.exe. Please add the folder containing it to PATH or set CUDA_ARCH_OVERRIDE in CMake.")
    	else()
    		set(NVIDIA_SMI_COMMAND "nvidia-smi")
            set(NVIDIA_SMI_FAILURE_MSG "nvidia-smi failed! Make sure the NVIDIA driver is installed and loaded!")
    	endif()
        execute_process(COMMAND ${NVIDIA_SMI_COMMAND} "-L"
                        RESULT_VARIABLE __res OUTPUT_VARIABLE __out
                        ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
    	if(__res EQUAL 0)
    		string(REGEX REPLACE "\n" ";" __out ${__out}) #Make into list
    		foreach(GPU_LINE ${__out})
    			string(REGEX REPLACE "GPU [0-9]: (.*) \\(UUID: .*" "\\1" GPU_NAME ${GPU_LINE})
    			string(REGEX REPLACE "GeForce " "" GPU_NAME ${GPU_NAME})
    			string(REGEX REPLACE "N[Vv][Ii][Dd][Ii][Aa] " "" GPU_NAME ${GPU_NAME})
    			get_compute_for_gpu(GPU_COMPUTE ${GPU_NAME})
    			message(STATUS "GPU found: " ${GPU_NAME} ", compute " ${GPU_COMPUTE})
    			set(TEMP "${TEMP} -gencode arch=compute_${GPU_COMPUTE},code=sm_${GPU_COMPUTE}")
    		endforeach()
    	else()
    		message(FATAL_ERROR "${NVIDIA_SMI_FAILURE_MSG}")
    	endif()
    endif()
    set(${OUT_VAR} ${TEMP} PARENT_SCOPE)
endfunction()

get_gencode_args(NVCC_GENCODE_FLAGS)
message(STATUS "Using gencode args: " ${NVCC_GENCODE_FLAGS})
set(CMAKE_CUDA_FLAGS "${NVCC_GENCODE_FLAGS} ${NVCC_DEBUG_FLAGS} ${CUDA_EXTRA_FLAGS}")
message(STATUS "Full CUDA flags: ${CMAKE_CUDA_FLAGS}")
