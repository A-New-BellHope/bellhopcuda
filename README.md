# bellhopcxx / bellhopcuda
C++/CUDA port of `BELLHOP`/`BELLHOP3D` underwater acoustics simulator.

### Impressum

Copyright (C) 2021-2022 The Regents of the University of California \
c/o Jules Jaffe team at SIO / UCSD, jjaffe@ucsd.edu \
Based on BELLHOP, which is Copyright (C) 1983-2022 Michael B. Porter

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <https://www.gnu.org/licenses/>.

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

- The performance of `bellhopcxx` roughly scales with the number of logical CPU
cores, as compared to `BELLHOP` which is single-threaded. On a 12-core /
24-thread test machine, `bellhopcxx` is typically about 10 to 30 times faster
than `BELLHOP`.
- For some very large simulations, `bellhopcuda` is even faster, occasionally 
approaching 100x gain over `BELLHOP`. For typical small simulations,
`bellhopcuda` is often slower than `bellhopcxx` multithreaded; see below.
- Our team has made [fixes to `BELLHOP`](https://github.com/A-New-BellHope/bellhop)
to improve its numerical properties and robustness (see below). These fixes were
essential to being able to reproduce its results in another language (or two).
These fixes are of course also incorporated into `bellhopcxx` / `bellhopcuda`.
- This repo also builds `bellhopcxx` as a library, which can be incorporated
into other programs, as well as the traditional executable version. The library
version allows the simulation parameters to be modified and the simulation
rerun, and it is not bottlenecked by slow file I/O for getting the results.

### Can I use `bellhopcxx` / `bellhopcuda` with MATLAB?

Yes. `bellhopcxx` and `bellhopcuda` read the same `.env` and other input files
and produce output files in the same formats as `BELLHOP`.

By default, one `bellhopcxx.exe` executable is produced, and you select the
dimensionality with command-line flags (`--2D`, `--3D`, `--Nx2D`, defaults to
2D). For no-code-changes MATLAB compatibility, you can set the `BHC_3D_SEPARATE`
CMake option, which will then instead produce separate `bellhopcxx2d.exe`,
`bellhopcxx3d.exe`, and `bellhopcxxnx2d.exe` executables, each only supporting
their respective dimensionality. (The corresponding CUDA executables are also
produced in each of these two modes, if enabled.)

Once you have your desired executable, simply rename it to `bellhop.exe` or
`bellhop3d.exe` and replace the `BELLHOP` / `BELLHOP3D` executable in your
MATLAB setup with it. Multithreading is enabled by default, so you should get a
significant speedup.

### Which features of the simulators are working?

Generally, all features of 2D and 3D are well-tested and reasonably mature (see
more detailed discussion below). In contrast, Nx2D has received little attention
so far; you can expect issues when trying to use Nx2D.

### How do I download the project?

We recommend `git clone` and build from source. You can [download pre-compiled
binaries for Windows 10 64-bit](https://github.com/A-New-BellHope/bellhopcuda/releases),
but these may be outdated compared to recent changes to the repo.

### How do I build the project?

- Make sure you got the Git submodules (the `glm` folder should not be empty).
- If you want `bellhopcuda`, install [the latest CUDA toolkit](https://developer.nvidia.com/cuda-downloads)
(see below for version information). Otherwise if you do not want `bellhopcuda`,
set the environment variable `BHC_NO_CUDA` (to something like "1") or turn off
the CMake option `BHC_ENABLE_CUDA`.
- Build the project with CMake in the usual way for your platform. If you are
not familiar with CMake, there are numerous tutorials online; the process is
not project-specific, besides the options mentioned in this readme.

### Why do I get errors about `hypot` when building with CUDA?

If you get an error like `error: ‘float cuda::std::__4::hypot(float, float, float)’
conflicts with a previous declaration` or similar, this is a bug in libcu++
which was fixed between libcu++ 1.7.0 (included in CUDA Toolkit 11.6) and
libcu++ 1.8.0. At time of writing, there is no CUDA toolkit download including
the latter, but if you are reading this in the future probably CUDA Toolkit
11.7 or later will include the fix. Until then, you have to 
[download libcu++ source manually](https://github.com/NVIDIA/libcudacxx)
and replace `include/cuda` and `include/nv` in your CUDA installation with the
appropriate directories from the Git repo.

### How do I use `bellhopcxxlib` as a library in another project?

- Copy the `[bellhopcuda repo]/include/bhc` directory to an `include` directory
in your other project. Or, if a Git submodule, add the `[bellhopcuda repo]/include`
directory to your include paths.
- Set your parent project's C++ version to C++17 or later.
- If you are on Windows and using the `bellhopcxxlib.dll` shared (dynamic) library
version, define `BHC_DLL_IMPORT` before including the header. This cannot be
done automatically because it must not be defined if you're using the static
library version, and we can't know in advance which version you will use.
- `#include <bhc/bhc.hpp>`. Follow the instructions in that file for API use.
- Link to:
    - Windows shared (dynamic) library: `bellhopcxxlib.dll` and import library `bellhopcxxlib.lib`
    - Windows static library: `bellhopcxxstatic.lib`
    - Linux shared library: `libbellhopcxxlib.so`
    - Linux static library: `libbellhopcxxstatic.a`

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

# Results

This section was last updated 5/2022; things may have changed since then.

## Accuracy

`bellhopcxx` / `bellhopcuda` includes a semi-automated test system. Tests can be
run automatically for all four run types (ray, transmission loss, eigenrays,
arrivals), and all results types are automatically checked, though the
capabilities of the automatic checkers vary for each type. The automatic
checkers check results based on absolute error, relative error, and ULPs, but
they don't use any contextual information. So for example a receiver being off
by 3\% at a level of 1e-7 is considered the same error whether the peak of the
field is 1e-7 around this receiver (this is a big deal), or the peak of the
field is 1e-2 hundreds of kilometers away (this is negligible).

The version of `BELLHOP` used when generating reference outputs to compare
results to should be
[our modified version of BELLHOP](https://github.com/A-New-BellHope/bellhop).
Many results of the original `BELLHOP` cannot be reproduced, in some cases not
even by `BELLHOP` itself; for more details, see the discussion on that repo.

Note that rays and eigenrays are generated in a random order due to thread
scheduling by `bellhopcxx` in multithreading mode and by `bellhopcuda`. Thus,
results comparison is sometimes only made to `bellhopcxx -1` (single-threaded
mode). However, note that the path of each ray considered independently does not
differ significantly between C++ and CUDA, as they are generated by the same
code. The only slight differences which may appear are due to floating-point
differences, e.g. CUDA merging floating-point multiplies and adds into fused
multiply-adds whenever possible.

#### Rays

The files `ray_tests_pass` and `ray_tests_fail` represent tests in which *both*
programs succeed or fail, respectively. That is, the latter are environment
files with formatting mistakes or other issues. Some of these environment files
may have been intended for other Acoustics Toolbox programs (not `BELLHOP`).

All the ray files produced in `ray_tests_pass` match, except for `MunkB_ray_rot`
which has one extra infinitesimal step on one ray in `BELLHOP` and a handful of
steps with around 0.0001% - 0.01% error.

#### Transmission Loss

`bellhopcxx` / `bellhopcuda` matches the results of `BELLHOP` on the 45 test
cases in `tl_match`, and fails to match on the three test cases in `tl_nomatch`.
In each of the latter cases, errors are restricted to one ray or one receiver,
and occur because of slight floating-point differences.

`aetB_TL` and `aetB` in `tl_long` produce results which have substantial errors,
but they take so long to run (about 41 minutes in `BELLHOP`) that they have not
been debugged yet.

#### Eigenrays

`bellhopcxx` / `bellhopcuda` matches the results of `BELLHOP` for the single
eigenrays test case provided by `BELLHOP`.

#### Arrivals

`bellhopcxx` / `bellhopcuda` matches the results of `BELLHOP` on the 6 test
cases in `arrivals_match` and fails to match on `sduct_bbB` in
`arrivals_nomatch`. In the latter, the problems are mild, affecting 19 out of
about 150,000 arrivals.

## Speed

`bellhopcxx` performance in multithreaded mode roughly scales with the number of
logical CPU cores. This is dependent on the problem size as well as other
bottlenecks such as file I/O for runs producing large output data files.

`bellhopcuda` performance for very large problem sizes does often outpace
`bellhopcxx` multithreaded. However, for most `BELLHOP` test files,
`bellhopcuda` is substantially slower than `bellhopcxx` multithreaded. This is
primarily for two reasons:
1. `bellhopcuda` has been tested on consumer ("gaming") GPUs only (e.g. NVIDIA
GeForce RTX 3080). These GPUs only include vestigial double-precision floating-
point support, at 1/64th the throughput compared to single-precision. Most of
`BELLHOP` and consequently `bellhopcxx` / `bellhopcuda` uses double-precision
(see below).
2. The codebase parallelizes over rays, the same for multithreaded C++ as for
CUDA. Unless there are many more rays than CUDA cores (at least several tens of
thousands for a high-end RTX GPU), the GPU will not be fully utilized.
Dedicated GPU-only code could be written with multiple threads contributing to
the same ray, e.g. by checking it against multiple receivers in parallel, which
could improve the performance.

`bellhopcxx` / `bellhopcuda` includes the `USE_FLOAT` option in CMake, which
forces the entire project to use single-precision (32-bit) floats.
This generally works but is not robust; for example, in a ray run, 99% of the
rays might follow the same trajectories as in the double version, just being off
by an infinitesimal amount; but 1% of the rays might diverge from their double
versions due to initially small numerical differences being amplified.

# Miscellaneous

## Comments

Unattributed comments in all translated code are copied directly from the
original `BELLHOP` and/or Acoustics Toolbox code, mostly by Michael B. Porter.
Unattributed comments in new code are by the Jules Jaffe team, mostly Louis
Pisha. It should usually be easy to distinguish the comment author from the
style.

## Code style

Initially, code style (things like where newlines and if statements are) was
kept as close as possible to the original Fortran code, for ease of comparing
the source files. As the project has evolved and templated 2D/Nx2D/3D was added,
the code style and organization has changed somewhat.
