# FAQs

### What is C++/CUDA?

This is a single codebase which can be built as multithreaded C++ code for your
CPU, or as CUDA code for your NVIDIA GPU. **You can use the CPU version 
(bellhopcxx) even if you don't have an NVIDIA GPU.**

### What platforms does this run on?

`bellhopcxx` is compatible with all platforms (Linux, Windows, Mac), and
`bellhopcuda` is compatible with all platforms which support CUDA (Linux and
Windows).

### Why should I use `bellhopcxx` / `bellhopcuda` instead of `BELLHOP`?

- While the performance gains from multithreading are highly dependent on the
  environment file, in most cases the performance scales roughly with the number
  of logical CPU cores. On a 12-core / 24-thread test machine, `bellhopcxx` is
  typically about 10 to 30 times faster than `BELLHOP`.
- For simulations with many (tens of thousands or more) rays and a reasonably
  good desktop GPU, `bellhopcuda` is typically even faster, sometimes surpassing
  100x gain over `BELLHOP`.
- Our team has made [fixes to `BELLHOP`](https://github.com/A-New-BellHope/bellhop)
to improve its numerical properties and robustness (see below). These fixes were
essential to being able to reproduce its results in another language (or two).
These fixes are of course also incorporated into `bellhopcxx` / `bellhopcuda`.
- As well as the traditional executable version, this repo also builds
  `bellhopcxx` as a library, which can be incorporated into other programs. The
  library version allows the simulation parameters to be modified and the
  simulation rerun, and it is not bottlenecked by slow file I/O for getting the
  results. The library version also allows running multiple independent
  simulations with different parameters at the same time (there are no global
  variables).
- By default, our version supports certain combinations of parameters which are
  not supported by the original version, such as 2D Gaussian ray-centered
  influence and 3D PCHIP. A compile-time option can be set to artificially limit
  the feature set to the same as `BELLHOP` / `BELLHOP3D` for testing.

### Can I use `bellhopcxx` / `bellhopcuda` with MATLAB?

Yes. `bellhopcxx` and `bellhopcuda` read the same `.env` and other input files
and produce output files in the same formats as `BELLHOP` / `BELLHOP3D`. It is
a drop-in replacement.

The normal `bellhopcxx.exe` and `bellhopcuda.exe` executables select their
dimensionality with command-line flags (e.g. `--2D`, `--3D`, `--Nx2D`). The
build system also produces executables which support only one dimensionality,
`bellhopcxx2d.exe`, `bellhopcxx3d.exe`, and so on. You can use these if you
don't want to change the MATLAB wrapper code to add the flags.

Once you have your desired executable, simply rename it to `bellhop.exe` or
`bellhop3d.exe` and replace the `BELLHOP` / `BELLHOP3D` executable in your
MATLAB setup with it. Multithreading is enabled by default, so you should get a
significant speedup.

### Which features of the simulators are working?

All features of all three dimensionalities (2D, 3D, Nx2D) are working and
reasonably well-tested. The numerical stability and robustness of all of
these versions is significantly improved from the original `BELLHOP` /
`BELLHOP3D` (see more discussion below), so in any application where the
original version was sufficient, our versions should be as well. With that said,
the numerical stability and reproducibility is better in 2D than in 3D or Nx2D,
and slightly better in 3D than Nx2D.

### How do I download the project?

[Download pre-compiled binaries for Windows 10 64-bit](
https://github.com/A-New-BellHope/bellhopcuda/releases). We try to keep these
up-to-date with most changes to the repo's main branch, though this is not
guaranteed. However, the pre-compiled versions:
- do not include all the executable and library versions possible to build
- are only built for one GPU version (SM86, RTX 30x0)
- are Windows versions only

If you need any other versions, please `git clone` and build from source.

### How do I build the project?

- Make sure you got the Git submodules (the `glm` folder should not be empty).
- If you want `bellhopcuda`, install [the latest CUDA toolkit](
  https://developer.nvidia.com/cuda-downloads). Otherwise if you do not want
  `bellhopcuda`, set the environment variable `BHC_NO_CUDA` (to something like
  "1") or turn off the CMake option `BHC_ENABLE_CUDA`.
- Build the project with CMake in the usual way for your platform. If you are
not familiar with CMake, there are numerous tutorials online; the process is
not project-specific, besides the options mentioned in this readme.

### Why do I get compiler errors when building with CUDA on Windows?

This project uses [libcudacxx](https://github.com/NVIDIA/libcudacxx) for
creating a shared codebase for CPU and GPU. This open-source library is
generally very good, but has occasional compatibility issues. Some errors you
might see might mention some of the following items:
- `cuda::std::__4::hypot`
- `_LInf`, `_LNan`, `_LSnan`
- `__ATOMIC_RELEASE`, `ATOMIC_INT_LOCK_FREE`
- `__type_traits/disjunction.h`, `error C2210: '_First': pack expansions cannot
  be used as arguments...`

We have [submitted an issue to NVIDIA](https://github.com/NVIDIA/libcudacxx/issues/354)
about some of these issues, and they are in the process of being fixed. We also
put workarounds in place in the `bellhopcxx` / `bellhopcuda` codebase. However,
if you run into any of these problems, always first upgrade to the latest CUDA
version. If that does not solve the problem, let us know and also try the
following:
- Upgrade from Visual Studio 2017, to 2019 or 2022.
- Download the libcudacxx source manually and replace `include/cuda` and
`include/nv` in your CUDA installation with the appropriate directories from
the Git repo.

### How do I use `bellhopcxxlib` as a library in another project?

- Copy the `[bellhopcuda repo]/include/bhc` directory to an `include` directory
in your parent project. Or, if a Git submodule, add the `[bellhopcuda repo]/include`
directory to your include paths.
- Set your parent project's C++ version to C++17 or later.
- If you are on Windows and using the `bellhopcxxlib.dll` shared (dynamic) library
version, define `BHC_DLL_IMPORT` before including the header. This cannot be
done automatically because it must not be defined if you're using the static
library version, and we can't know in advance which version you will use.
- `#include <bhc/bhc.hpp>`. Follow the instructions in that file for API use.
- Link to:
    - Windows shared (dynamic) library: `bellhopcxxlib.dll` and import library
      `bellhopcxxlib.lib`
    - Windows static library: `bellhopcxxstatic.lib`
    - Linux shared library: `libbellhopcxxlib.so`
    - Linux static library: `libbellhopcxxstatic.a`

### What additional functionality is available in library mode?

- You can set up the library without an initial environment file. All parameters
are set to reasonable basic defaults. You can modify these up one at a time from
your host program. See examples in `examples`.
- The `extsetup` functions allow you to allocate arrays with the run parameters,
for example to change the number of sources in Z to 100. After allocating the
array with your chosen size, you can fill in the data directly. 
- The `writeenv` function allows you to save the current state of the run
parameters to an environment file and other required input files (e.g. SSP,
bathymetry). This does not have to be relative to the same FileRoot provided at
setup time.
- The `readout` function allows you to read results from a past run (ray file,
TL / shade file, or arrivals) into memory, so your host program can display or
manipulate these results.

### How do I report bugs?

If the bug is regarding results being wrong (or different from
`BELLHOP`/`BELLHOP3D`), please check and report to us whether the behavior is
the same in [our version of `BELLHOP`/`BELLHOP3D` with improved numerical
stability and robustness](https://github.com/A-New-BellHope/bellhop). Also,
please provide an environment file which traces only the ray whose results are
incorrect (if applicable); if this is impossible, due to the bug disappearing
when a different set of rays are traced, let us know of this.

Submit a bug report on the
[GitHub issues page](https://github.com/A-New-BellHope/bellhopcuda/issues).
