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

### Building bellhopcuda (with HIP/ROCm support, AMD GPUs)
To build the HIP/ROCm version for AMD GPUs, ensure you have
[ROCm installed](https://rocm.docs.amd.com/) and available in your system PATH.
Follow the CPU steps above but configure with `-DBHC_ENABLE_HIP=ON -DBHC_ENABLE_CUDA=OFF`
(HIP and CUDA are mutually exclusive). You may need to add
`-DCMAKE_HIP_COMPILER=/opt/rocm/llvm/bin/clang++` so CMake finds the HIP compiler.

The target GPU architecture is set with `-DCMAKE_HIP_ARCHITECTURES=<arch>` (e.g.
`gfx90a`, `gfx1100`), defaulting to `gfx90a` when unset. We have tested on AMD
Instinct (gfx90a, CDNA2) and Radeon (gfx1100 RDNA3, gfx1201 RDNA4) GPUs.

These compilation paths will produce a set of executables and libraries in a bin directory, with the
CUDA- and HIP-enabled versions (bellhopcuda*) having additional GPU support. Note that building with CUDA
is hardware specific; ensure your GPU is compatible with the CUDA version you have installed.
