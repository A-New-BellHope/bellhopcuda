# bellhopcxx / bellhopcuda
C++/CUDA port of `BELLHOP`/`BELLHOP3D` underwater acoustics simulator.

### Impressum

Copyright (C) 2021-2026 The Regents of the University of California \
Marine Physical Lab at Scripps Oceanography, c/o Jules Jaffe, jjaffe@ucsd.edu \
Based on BELLHOP / BELLHOP3D, which is Copyright (C) 1983-2022 Michael B. Porter

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <https://www.gnu.org/licenses/>.

## Breaking changes

While we strive to maintain backward compatibility with bellhopcxx /
bellhopcuda version, some breaking changes have been introduced in version 1.5.0.
Most notably, the binarry output data files, **.shd and .arr**, have been updated
to match the current behavior of [`BELLHOP` / `BELLHOP3D`](http://oalib.hlsresearch.com/AcousticsToolbox/)
2024 update to the Acoustics Toolbox.

CMake has also been updated to version 3.27 minimum to support automated recognition
of CUDA enabled systems. You may need to update your CMake installation.

## Updated Features

- Provence ocean SSP model support.
- Background calculation in library mode.
- Combined Eigenray and arrivals run mode using `V` option.
- Several bug fixes and performance improvements.
- More example programs showing [library usage](examples/).

## Current Status

There is an ongoing issue with 3D runs in complex environments falsely
identifying rays as interacting with receivers, even at very large distances.
This can lead to errors in all run types except ray runs. We are actively working
to resolve this issue. This issue affects both CPU
multithreaded and CUDA modes. It is also present in the original `BELLHOP` /
`BELLHOP3D` code.

## Compilation

We use a CMake-based build system. Please see the
[doc/compilation.md](doc/compilation.md) document for detailed instructions.

## Speed

- The speedup in multithreaded mode roughly scales with the number of
  logical CPU cores.
- The speedup in CUDA mode
  - Consumer-grade GPU (e.g. NVIDIA GeForce RTX 3060): typically about 10x-50x
	speedup over `BELLHOP` / `BELLHOP3D` for large runs with few receivers.
  - Server-grade GPU (e.g. NVIDIA A100): typically about 20x-100x speedup and 
	less sensitive to receivers layout.

Detailed performance results are available in the [doc/performance.md](doc/performance.md)
document.

## Accuracy

The accuracy characteristics of `bellhopcxx` / `bellhopcuda` are discussed in
the [doc/accuracy.md](doc/accuracy.md) document.

# Miscellaneous

## Comments

Unattributed comments in all translated code are copied directly from the
original `BELLHOP` / `BELLHOP3D` and/or Acoustics Toolbox code, mostly by Dr.
Michael B. Porter. Unattributed comments in new code are by the Marine Physical
Lab team, mostly Louis Pisha. It should usually be easy to distinguish the
comment author from the style.
