# bellhopcxx / bellhopcuda
C++/CUDA port of `BELLHOP`/`BELLHOP3D` underwater acoustics simulator.

### Impressum

Copyright (C) 2021-2023 The Regents of the University of California \
Marine Physical Lab at Scripps Oceanography, c/o Jules Jaffe, jjaffe@ucsd.edu \
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
and produce output files in the same formats as `BELLHOP` / `BELLHOP3D`.

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
- Upgrade from Visual Studio 2017 to 2019 or 2022.
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

This section was last updated 2/2023.

## Speed

The summary is:
- Speedups vary widely depending on many parameters of the run.
- Generally, the speedup in multithreaded mode roughly scales with the number of
  logical CPU cores, with a scaling factor roughly ranging from 0.5 to 1.2
  (e.g. a 6-core, 12-thread CPU will often get a 6x-15x speedup).
- If you have a reasonably modern (but consumer-grade) GPU and the run uses a
  large number of rays (in the 10,000s or more), the speedup of the CUDA version
  will often be a few times above the CPU speedup (e.g. 20x-100x compared to the
  Fortran version). Conversely, if there are a small number of rays, the CUDA
  performance will be worse than the CPU performance.

Some factors affecting the performance are discussed below.

#### Ray count

In both CPU multithreaded mode and CUDA mode, the processing is parallelized
over the rays. The number of rays needed to fully utilize each chip is a few
times the number of cores--but on the CPU this will be anything over a hundred
or so rays, whereas on the GPU you will need tens of thousands of rays as there
are many thousands of cores. However, experimentally, the speedup does not stop
once this number is reached--running hundreds of thousands or millions of rays
continues to bring additional, though diminishing, performance gain over
`BELLHOP` / `BELLHOP3D`. Of course, whether more rays is useful or not is very
application- dependent.

#### File I/O

`bellhopcxx` / `bellhopcuda` can be built and used as a library, so input data
(e.g. SSP) does not have to be read from files, and output data does not have to
be written to disk. In some runs, the file I/O is trivial, whereas in others,
it takes 10x-100x the time compared to actually computing the results. Normally,
we count the file I/O time in our version in order to compare most fairly to
the original version, but if file I/O can be avoided in a particular
application using the library version, this may have a substantial impact on the
performance.

#### Receivers layout

A run with receivers concentrated in one area far away from the source will be
faster than one where the receivers are spread uniformly over a large area where
the rays traverse. In the latter, every step of every ray will influence a
number of receivers, whereas in the former, whereas in the former most rays will
not influence any receivers. This is especially true for eigenrays and arrivals
runs: these run types will work best when the user is interested in finding a
handful of rays which reach an area of interest out of a large number of rays
initially traced.

#### Run type

Transmission loss runs typically bring the largest speedups over `BELLHOP` /
`BELLHOP3D`. Arrivals typically have lower speedups, depending on the receivers'
layout as discussed above.

For ray runs, the full rays must be written out to disk, serially one after
another. For this reason, our CUDA version does not even compute the rays on the
GPU at all: it is unlikely the user will want to run tens of thousands of rays
to saturate the GPU, and if they did, writing them is likely to take more time
than computing them. Nevertheless, the performance gains with multithreaded mode
are still often noticeable.

For eigenrays runs, when a ray influences a receiver, a data structure similar
to an arrival is written to memory, but the full ray is not saved to memory (or
worse for the parallelism, to disk). This allows a large number of rays to be
searched for the few which hit target receivers, on either CPU multithreaded or
CUDA. Then, at the end, just those rays which hit receivers are re-traced in CPU
multithreaded mode (even if the CUDA version is running), and the resulting rays
are written. This can be much more efficient than `BELLHOP` / `BELLHOP3D` when
there are few receivers, but it will be less efficient when there are receivers
everywhere and consequently extremely large numbers of eigenrays.

#### Dimensionality

The speedup for 2D is typically slightly (maybe 10% on average) higher than for
3D. The speedup for Nx2D is similar to that for 3D.

#### Floating point precision

By default, `bellhopcxx` / `bellhopcuda` uses double precision (64-bit floats)
for almost all arithmetic. (It uses the same precision as the Fortran version
for each individual operation, which is almost universally double precision,
except for most writing of output data which is in single precision.) All the
performance and accuracy results presented here are for double precision mode.

However, our version can also be built to use single precision (32-bit floats)
for all arithmetic. This generally works--the ray or TL plots generally look the
same by eye--but as you might expect, this substantially increases the chances
of numerical stability problems, and the details of the results generally do not
match the double precision version. For example, in a particular run, 99% of the
rays might follow nearly the same trajectory as with full precision, but 1% of
the rays might diverge and go elsewhere.

Consumer NVIDIA GPUs only contain vestigial double-precision floating-point
support, at 1/64th the throughput compared to single-precision. Only the server-
class GPUs (in the tens of thousands of dollars) include full double-precision
support. The fact that the CUDA version still can bring substantial speedups
over the CPU is explained by the fact that most of the actual instructions being
executed are integer or memory instructions, with only a portion being
floating-point math. Nevertheless, the floating-point performance may be the
bottleneck in some runs. For some applications, the single-precision version
may be useful for obtaining fast, approximate initial results on a consumer GPU. 

## Accuracy

The physics model in the original `BELLHOP` / `BELLHOP3D` has a number of
properties which make it sensitive to numerical details. For example, moving a
source by 1 cm in a basic test case causes up to about 3 dB differences
throughout the field. As a more extreme example, we created a modified version
of the original `BELLHOP` code which adds random perturbations less than *100
microns* to the position (not even direction!) of each ray after each step.
(Steps are typically on the order of 1 km, so this is a relative error on the
order of 10^-7.) This was chosen because on the one hand, the real-world
uncertainty in the ocean is at a far larger scale, and on the other hand, these
perturbations wreak havoc on `BELLHOP`'s edge case handling (see below). These
tiny perturbations result in up to 40 dB differences in the field, which are
visible when comparing the transmission loss plots by eye.

Of course, `BELLHOP` has been providing useful results to researchers for four
decades now, despite these limitations of its model. We mention this for the
sake of context for the discussion below. Usually in software development, any
mismatch between the program's output and the reference output means the test
has failed and the program has a bug somewhere. In this project, such a mismatch
could be caused by something as subtle as the Fortran compiler emitting a fused
multiply-add instruction where the C++ compiler emitted separate multiply and
add instructions. (In fact, in one case, we confirmed this was the source of a
discrepancy from examining the disassembly, and explicitly added the FMA to the
C++ version.)

But not only is it not realistic to try to reproduce the exact set of floating-
point operations in a different language and compiler to get exactly matching
results--we *want* the set of operations to be different in some cases. For
example, early on in the project, we found a combination of bugs and edge case
conditions--documented in [the readme of our modified Fortran
version](https://github.com/A-New-BellHope/bellhop)--which made it impossible to
parallelize the computation if we wanted to reproduce the original results (with
the bugs). Instead of abandoning the parallelization--which was one of the main
goals of the project--we decided to fix the bugs and improve the numerical
stability.

Since then, we have been comparing the results of `bellhopcxx` / `bellhopcuda`
to our modified Fortran version, not to the original `BELLHOP` / `BELLHOP3D`.
When we find a discrepancy, we find the failing ray, print out its state as it
travels, find the step where it diverges, and isolate the part of the code
responsible. Of course, if it is a bug in our new code, we fix that. But when it
is a numerical stability / robustness issue, we implement a change in *both*
versions which attempts to make that part of the code more reproducible, and
hopefully then the results start matching again.

We have made a variety of these changes, [documented
here](https://github.com/A-New-BellHope/bellhop), but the most significant set
of changes was to the handling of boundaries. Boundaries are locations in the
ocean, such as a change in SSP or a change in the bathymetry slope, where the
ray may have to recompute its direction to curve. As such, the ray steps forward
until it hits the next boundary in the direction it is traveling, and then
adjusts its direction. Thus, most steps place a ray directly on a boundary--in
other words, almost every step of every ray hits an edge case. Due to
floating-point imprecision, this may actually be slightly in front of or behind
the boundary. This affects the subsequent trajectory in a way which is usually
very small, but can be amplified to large changes later. The original `BELLHOP`
/ `BELLHOP3D` made no attempt to handle stepping to boundaries consistently; our
version handles most boundaries in a fully consistent manner and the remaining
couple in a manner which is perhaps about 99.99% consistent. (Of course, that
means the remaining ~0.01% of rays diverge slightly, potentially causing
results to not match.)

So all of the results discussed here are as compared to [our modified Fortran
version](https://github.com/A-New-BellHope/bellhop), not the original `BELLHOP`
/ `BELLHOP3D`. And when we report results that do not match, these are mostly
either cases where the best methods we were able to come up with still
occasionally diverge, or cases we have not studied in detail. Because of how
many cases do match, it's unlikely that these are just the result of bugs in our
version, though it is possible that there are bugs on codepaths which have been
rarely or never tested such as some of the more obscure attenuation options.

One final note: We generally compare values, whether ray positions or complex
field values, with a combination of absolute and relative error thresholds.
There is also logic for things such as reordering and getting rid of duplicate
steps in rays, and combining arrivals. However, the criteria are always much
more strict than what would be noticed from simply looking at MATLAB plots of
the outputs.

### Coverage Tests

A set of scripts are provided which generate test environment files (along with
SSP files, bathymetry, etc. as appropriate) for every combination of run type,
influence type, SSP, and dimensionality (about 1500 tests). These are primarily
for code coverage, rather than testing the correctness of the physics, and are
short to be able to run relatively quickly. The scripts also test to make sure
that any combinations of these settings which are not supported by `BELLHOP` /
`BELLHOP3D` are rejected by our version, assuming the `BHC_LIMIT_FEATURES` build
option was enabled.

The current version of `bellhopcxx` / `bellhopcuda` matches on all the coverage
tests.

### Dr. Michael B. Porter's Tests

Dr. Porter provided a set of test cases with the Acoustics Toolbox. These have
been sorted by dimensionality and run type. The specific runs comprising this
data are in the text files included in the root of this repo, e.g.
`ray_3d_pass.txt`, `tl_long.txt`, etc. This table is just counts of the entries
in those files.

The "Both fail" column is for tests which fail to run in `BELLHOP` / `BELLHOP3D`
at all, either because they may be for one of the other Acoustics Toolbox
simulators or because the environment file syntax has changed since they were
written. Our version produces the same behavior as the original version on all
of these test cases (though often with more descriptive error messages), so this
"Both fail" case is actually a "pass" for our version.

|     | Match | Do not match | Both fail |
| --- | --- | --- | --- |
| 2D ray | 12 | 0 | 8 |
| 2D TL | 48 | [1] | 0 |
| 2D eigen | 1 | 0 | 0 |
| 2D arrivals | 6 | 1 | 0 |
| 3D ray | 12 | 0 | 2 |
| 3D TL | 42 | 8 | 12 |
| 3D eigen | 0 | 0 | [2] |
| 3D arrivals | 0 | 0 | 0 |
| Nx2D ray | 4 | 0 | 0 |
| Nx2D TL | 13 | 5 | 6 |
| Nx2D eigen | 0 | 0 | 0 |
| Nx2D arrivals | 0 | 0 | 0 |

[1]: There are two runs which consistently match in single-threaded mode, but
sometimes do and sometimes don't match in multithreaded or CUDA mode. \
[2]: There is only one environment file; it runs out of memory after over 3
hours, so this has not been investigated further.

Note that some of the TL tests take several hours (one about 12 hours), and we
do not guarantee we will re-run those tests after every change, so some of the
numbers here may be slightly out of date.

# Miscellaneous

## Comments

Unattributed comments in all translated code are copied directly from the
original `BELLHOP` and/or Acoustics Toolbox code, mostly by Michael B. Porter.
Unattributed comments in new code are by the Jules Jaffe team, mostly Louis
Pisha. It should usually be easy to distinguish the comment author from the
style.
