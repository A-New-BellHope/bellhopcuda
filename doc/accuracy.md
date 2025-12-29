# Accuracy of bellhopcxx / bellhopcuda

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
goals of the project--we decided to fix (many of) the bugs and improve the
numerical stability.

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
