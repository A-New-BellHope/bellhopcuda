# Compilation of bellhopcxx / bellhopcuda

We use a CMake-based build system. The following instructions assume you have
[CMake >=3.27 installed](https://cmake.org/install/) and available in your system PATH.
Build compiles and runs on Linux and Windows, and we recommend the most recent compilers with
C++17 support.

### Building bellhopcxx (CPU only)
To build the CPU-only version, `bellhopcxx`, follow these steps:

1. Clone the repository:
   ```bash
   git clone https://github.com/A-New-BellHope/bellhopcuda.git
   cd bellhopcuda
   git submodule update --init
   git submodule update --recursive
   ```
2. Create a build directory and navigate into it:
   ```bash
   mkdir build
   cd build
   cmake .. -DBHC_ENABLE_CUDA=OFF
   cmake --build .
   ``` 

On Windows using Visual Studio 2022 or later, you can open the folder in Visual Studio directly after
cloning, and it will automatically configure the project.

### Building bellhopcuda (with CUDA support)
To build the CUDA-enabled version, `bellhopcuda`, ensure you have the
[NVIDIA CUDA Toolkit installed](https://developer.nvidia.com/cuda-downloads) and
available in your system PATH. Then follow the steps above , but omit the
`-DBHC_ENABLE_CUDA=OFF` option in the `cmake` command.

We have tested on many GPUs, including consumer models from the 20xx, 30xx, and 40xx series, server
GPUs A6000, A100, and GH200. We recommend using the latext version of CUDA and commonly compile
on CUDA versions 12.4, and 12.8, and 13.1.

Both compilation paths will produce a set of executables and libraries in a bin directory, with the
CUDA-enabled version (bellhopcuda*) having additional GPU support. Note that building with CUDA
is hardware specific; ensure your GPU is compatible with the CUDA version you have installed.
