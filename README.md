# bellhopcuda
CUDA and C++ port of BELLHOP.

## Comments
Unattributed comments in all translated code are copied directly from the original
BELLHOP and/or Acoustics Toolbox code, mostly by Michael B. Porter. Unattributed
comments in new code are by the bellhopcuda team, mostly Louis Pisha.

## Code style
Code style (things like where newlines and if statements are) is kept as close
as possible to the original FORTRAN code, for ease of comparing the source files.

# Map

### Support code

`common.hpp`: Main support code, CUDA compatibility, utilities, etc.
- `misc/monotonicMod.f90`
- `misc/SortMod.f90`
- `misc/MathConstants.f90`
- `misc/cross_products.f90`: Only used by BELLHOP3D, provided by glm.

`ldio.hpp`: C++ emulation of FORTRAN list-directed I/O.

`subtab.hpp`: Feature that allows a list to be automatically linearly interpolated
between its extreme values.
- `misc/subtabulate.f90`

`curves.hpp`: Templated splines and Piecewise Cubic Hermite Interpolating Polynomial.
- `misc/pchipMod.f90`: NEED all
- `misc/splinec.f90`: NEED all (templated version for real/cpx)
- `misc/splined.f90`
- `misc/spliner.f90`


### Data I/O code

`ssp.hpp`: Sound speed as a function of position.
- `Bellhop/sspMod.f90`
- `misc/sspMod.f90`: Base implementation for other programs, but without
derivatives and with less pre-computed info. Not used by BELLHOP.

`attenuation.hpp`: Sound speed and attenuation.
- `misc/AttenMod.f90`

`boundary.hpp`: Ocean surface and floor.
- `Bellhop/bdryMod.f90`: NEED all
- `Bellhop/bdry3DMod.f90`: BELLHOP3D version.

`refcoef.hpp`: Reflection coefficients.
- `misc/RefCoef.f90`: InterpolateReflectionCoefficient, NEED ReadReflectionCoefficient

`sourcereceiver.hpp`: Source and receiver (single or array) positions.
- `misc/SourceReceiverPositions.f90`

`beampattern.hpp`: Source beam pattern.
- `misc/beampattern.f90`: NEED all

`readenv.cpp`: Main environment file reading.
- `Bellhop/ReadEnvironmentBell.f90`, NEED OpenOutputFiles
- `misc/ReadEnvironmentMod.f90`: Base implementation for other programs, not
used by BELLHOP.
- `misc/RWSHDFile.f90`: NEED all

### Core simulation code

`structs.hpp`: Misc. BELLHOP structs not in other files.

`step.hpp`: Ray/beam tracing single step.
- `Bellhop/Step.f90`
- `Bellhop/Step3DMod.f90`: BELLHOP3D version.

`trace.hpp`: Ray/beam tracing main.
- `Bellhop/bellhop.f90`: TraceRay2D, Distances2D, Reflect2D
- `Bellhop/bellhopMod.f90`: structs
- `Bellhop/ReflectMod.f90`: Alternate implementation of Reflect2D. BELLHOP
seems to use the implementation within `Bellhop/bellhop.f90`. BELLHOP3D uses
this (despite it not being 3D).
- `Bellhop/Reflect3DMod.f90`: BELLHOP3D version.

### Top-level executables

`Bellhop/bellhop.f90`: 
`Bellhop/bellhop3D.f90`: 








### FORTRAN code without C++ equivalents:
- `misc/FatalError.f90`: Not applicable.

### FORTRAN code not used by BELLHOP or BELLHOP3D:
- `misc/calculateweights.f90`
- `misc/interpolation.f90`
- `misc/MergeVectorsMod.f90`
- `misc/munk.f90`
- `misc/norms.f90`
- `misc/PekRoot.f90`
- `misc/PolyMod.f90`
- `misc/RootFinderSecantMod.f90`

`Bellhop/angleMod.f90`: 
`Bellhop/ArrMod.f90`: 


`Bellhop/Cone.f90`: 
`Bellhop/influence.f90`: 
`Bellhop/influence3D.f90`: 
`Bellhop/RayNormals.f90`: 
`Bellhop/WriteRay.f90`: 

`misc/weight New in progress.f90`
