# Performance of bellhopcxx / bellhopcuda

Some factors affecting the performance are discussed below.

## Ray count

In both CPU multithreaded mode and CUDA mode, the processing is parallelized
over the rays. The number of rays needed to fully utilize each chip is a few
times the number of cores--but on the CPU this will be anything over a hundred
or so rays, whereas on the GPU you will need tens of thousands of rays as there
are many thousands of cores. However, experimentally, the speedup does not stop
once this number is reached--running hundreds of thousands or millions of rays
continues to bring additional, though diminishing, performance gain over
`BELLHOP` / `BELLHOP3D`. Of course, whether more rays is useful or not is very
application-dependent.

## File I/O

`bellhopcxx` / `bellhopcuda` can be built and used as a library, so input data
(e.g. SSP) does not have to be read from files, and output data does not have to
be written to disk. In some runs, the file I/O is trivial, whereas in others,
it takes 10x-100x the time compared to actually computing the results. Normally,
we count the file I/O time in our version in order to compare most fairly to
the original version, but if file I/O can be avoided in a particular
application using the library version, this may lead to a substantial gain in
effective overall performance in some cases.

## Receivers layout

A run with receivers concentrated in one area far away from the source will be
faster than one where the receivers are spread uniformly over a large area where
the rays traverse. In the latter, every step of every ray will influence a
number of receivers, whereas in the former, whereas in the former most rays will
not influence any receivers. This is especially true for eigenrays and arrivals
runs: these run types will work best when the user is interested in finding a
handful of rays which reach an area of interest out of a large number of rays
initially traced.

![Receiver layout effect on performance](images/receiver_layout_performance.png)]

The receiver layout limits the performance of the GPU model, especially on consumer
GPUs with less memory bandwidth. Each ray must check at each step whether it influences
any receivers, and if so, read and write data to memory for each influenced receiver. On the
GPU, many rays are being processed simultaneously, and if many of them influence
receivers at the same time, this can lead to memory access bottlenecks which
limit the performance. On server-grade GPUs with higher memory bandwidth, this
effect is less pronounced.

## Run type

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
