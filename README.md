# bellhopcxx / bellhopcuda
C++/CUDA port of `BELLHOP` underwater acoustics simulator

### Impressum

Copyright (C) 2021-2022 The Regents of the University of California
c/o Jules Jaffe team at SIO / UCSD, jjaffe@ucsd.edu
Based on BELLHOP, which is Copyright (C) 1983-2020 Michael B. Porter

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
CPU, or as CUDA code for your NVIDIA GPU.

### Why should I use `bellhopcxx` / `bellhopcuda` instead of `BELLHOP`?

- The performance of `bellhopcxx` roughly scales with the number of logical CPU
cores, as compared to `BELLHOP` which is single-threaded. On a 12-core /
24-thread test machine, `bellhopcxx` is typically about 10 to 30 times faster
than `BELLHOP`.
- For some very large simulations, `bellhopcuda` is even faster, occasionally 
approaching 100x gain over `BELLHOP`. (For typical small simulations,
`bellhopcuda` is often slower than `bellhopcxx` multithreaded; see below.)
- Our team has made [fixes to `BELLHOP`](https://github.com/A-New-BellHope/bellhop)
to improve its numerical properties and robustness (see below). These fixes were
essential to being able to reproduce its results in another language (or two).
These fixes are of course also incorporated into `bellhopcxx` / `bellhopcuda`.
- This repo also builds `bellhopcxx` as a library, which can be incorporated
into other programs, as well as the traditional executable version. The library
version allows the simulation parameters to be modified and the simulation
rerun, and it is not bottlenecked by slow file I/O for getting the results.

### Can I use `bellhopcxx` / `bellhopcuda` with MATLAB?

Yes, simply rename `bellhopcxx.exe` to `bellhop.exe` and replace the Fortran
`BELLHOP` executable in your MATLAB setup with it. It reads the same `.env` and
other input files and produces output files in the same formats as `BELLHOP`.

### How do I download the project?

We recommend `git clone` and build from source. You can [download pre-compiled
binaries for Windows 10 64-bit](https://github.com/A-New-BellHope/bellhopcuda/releases),
but these may be outdated compared to recent changes to the repo.

### How do I build the project?

- Make sure you got the Git submodules (the `glm` folder should not be empty).
- If you want `bellhopcuda`, install [the latest CUDA toolkit](https://developer.nvidia.com/cuda-downloads)
(the earliest supported CUDA version is somewhere between 11.3 and 11.5).
Otherwise if you do not want `bellhopcuda`, set the environment variable
`BHC_NO_CUDA` or turn off the CMake option `BHC_ENABLE_CUDA`. 
- Build the project with CMake as usual.

### How do I use `bellhopcxxlib` as a library in another project?

- Copy the `[bellhopcuda repo]/include/bhc` directory to an `include` directory
in your other project. Or, if a Git submodule, add the `[bellhopcuda repo]/include`
directory to your include paths.
- `#include <bhc/bhc.hpp>`. Follow the instructions in that file for API use.
- Link to `bellhopcxxlib.dll` / `libbellhopcxxlib.so`.

### How do I report bugs?

If the bug is regarding results being wrong (or different from `BELLHOP`),
please check and report to us whether the behavior is the same in [our version
of `BELLHOP` with improved numerical stability and robustness](https://github.com/A-New-BellHope/bellhop).
Also, please provide an environment file which traces only the ray whose results
are incorrect (if applicable); if this is impossible, due to the bug
disappearing when a different set of rays are traced, let us know of this.

Submit a bug report on the [GitHub issues page](https://github.com/A-New-BellHope/bellhopcuda/issues).

# Results

This section updated 3/2022.

## Accuracy

`bellhopcxx` / `bellhopcuda` includes a semi-automated test system. Tests are
run automatically for all four run types (ray, transmission loss, eigenrays,
arrivals), but only transmission loss and arrivals results are automatically
checked. The automatic checkers check results based on absolute error, relative
error, and ULPs, but they don't use any contextual information, so for example
a receiver being off by 1\% at a level of 1e-7 is considered the same error
whether the peak of the field is 1e-7 around this receiver, or the peak of the
field is 1e-2 hundreds of kilometers away.

All results are compared to the outputs of [our modified version of BELLHOP](https://github.com/A-New-BellHope/bellhop).
Many results of the original BELLHOP cannot be reproduced, some not even by
itself; for more details, see the discussion on that repo.

Note that rays and eigenrays are generated in a random order due to thread
scheduling by `bellhopcxx` in multithreading mode and by `bellhopcuda`. Results
comparison is made to `bellhopcxx -1` (single-threaded mode). However, the
path of each ray considered independently does not differ significantly between
C++ and CUDA, as they are generated by the same code. The only slight
differences which may appear are due to floating-point differences, e.g. CUDA
merging floating-point multiplies and adds into fused multiply-adds whenever
possible.

#### Rays

All the ray tests provided in `BELLHOP` match in `bellhopcxx` / `bellhopcuda` as
evaluated by eye.  The files `ray_tests_pass` and `ray_tests_fail` represent
tests in which *both* programs succeed or fail, respectively. (That is, the
latter are environment files with formatting mistakes or other issues.)

All the ray files produced in `ray_tests_pass` have the same number of steps,
except that `BELLHOP` produces an extra infinitesimal step at line 866 of
`MunkB_ray_rot`. Values are typically accurate to at least 7-8 significant
figures.

#### Transmission Loss

`bellhopcxx` / `bellhopcuda` matches the results of `BELLHOP` on the 45 test
cases in `tl_match`, and fails to match on the three test cases in `tl_nomatch`.
In each of the latter cases, errors are restricted to one ray or one receiver,
and occur because of slight floating-point differences.

`aetB_TL` and `aetB` in `tl_long` produce results which are completely wrong,
but they take so long to run (about 41 minutes in `BELLHOP`) that they have not
been debugged yet.

#### Eigenrays

`bellhopcxx` / `bellhopcuda` matches the results of `BELLHOP` for the single
eigenrays test case provided by `BELLHOP`.

#### Arrivals

`bellhopcxx` / `bellhopcuda` matches the results of `BELLHOP` on the 6 test
cases in `arrivals_pass` and fails to match on `sduct_bbB` in `arrivals_fail`.
In the latter, the problems are mild, affecting 19 out of about 150,000
arrivals.

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
`BELLHOP` and consequently `bellhopcuda` uses double-precision (see below).
2. The codebase parallelizes over rays, the same for multithreaded C++ as for
CUDA. Unless there are many more rays than CUDA cores (at least several tens of
thousands for a high-end RTX GPU), the GPU will not be fully utilized.
Dedicated GPU-only code could be written with multiple threads contributing to
the same ray, e.g. by checking it against multiple receivers in parallel, which
could improve the performance.

`bellhopcxx` / `bellhopcuda` includes the `USE_FLOAT` option in CMake, which
forces the entire project to use single-precision (32-bit) floats.
Unfortunately, this leads to massive numerical precision issues, and the results
are unusable. Improving the algorithms' numerical properties to permit
generation of reasonable results in single-precision mode is one possible future
research direction.

# Miscellaneous

## Comments
Unattributed comments in all translated code are copied directly from the original
BELLHOP and/or Acoustics Toolbox code, mostly by Michael B. Porter. Unattributed
comments in new code are by the Jules Jaffe team, mostly Louis Pisha. It should
usually be easy to distinguish the comment author from the style.

## Code style
Code style (things like where newlines and if statements are) is kept as close
as possible to the original Fortran code, for ease of comparing the source files.
