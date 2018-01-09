# FourierFlows.jl Documentation


## Overview

FourierFlows provides solvers for partial differential equations on
doubly-periodic domains using Fourier-based pseudospectral methods.
A central intent of the software's design is also to provide a framework
for writing new, fast solvers for new physical problems.
The code is written in [Julia](https://julialang.org).


## Installation

FourierFlows is a registered package and so can be installed via `Pkg.add`.

```julia
Pkg.add("FourierFlows")
```

For now, this package supports Julia `0.6`.


## Source code organization

The code is divided along conceptual lines into problem-agnostic and
problem-specific components. Files that contain problem-agnostic parts
of the code are stored in `/src`. Files in `/src` define the domain,
'AbstractTypes' that supertype problem-specific types, and
time-stepper types and routines. Problem-specific modules are stores
in `/src/physics`.

Here's an overview of the code structure:

- `/src/`
    - `FourierFlows.jl`
        - Defines supertyping AbstractParams, AbstractGrid, etc.
        - Defines a `Problem` type to organize the grid, vars, params,
            equation, and timestepper into a single structure.
        - Includes all sources files and physics files.
   - `timesteppers.jl`: defines modules and `stepforward!` routines for
        various time-steppers. Current implemented time-steppers are:
        - Forward Euler (+ Filtered Forward Euler)
        - 3rd-order Adams-Bashforth (AB3)
        - 4th-order Runge-Kutta (RK4)
        - 4th-order Runge-Kutta Exponential Time Differencing (ETDRK4)
        (+ Filtered ETDRK4)
    - `physics/`
        - `twodturb.jl`: Defines a `TwoDTurb` module that provides a
                solver for the two-dimensional vorticity equation.
        - `barotropicqg.jl`: Defines a `BarotropicQG` module that provides
                several solvers for the barotropic QG model that permit beta,
                topography, beta + topography, and forcing.
        - `twomodeboussinesq.jl`: Defines a `TwoModeBoussinesq` module
                that provides solvers for a two-mode truncation of the
                rotating, stratified Boussinesq equation.
        - `niwqg.jl`: Defines a `NIWQG` module that provides a solver
                for the vertical-plane-wave model for the interaction of
                a near-inertial wave field and quasi-geostrophic flow.
        - `traceradvdiff.jl`: Defines a `TracerAdvDiff` module that
                provides a solver for a two-dimensional and periodic tracer
                field in a given 2D flow (u, w), which can be an arbitrary
                function of x, z, and t.


## Writing fast solvers

The performance-intensive part of the code involves just two functions: the
timestepping scheme `stepforward!`, and the function `calcNL!` that
calculates the nonlinear part of the given equation's right-hand side.
Optimization of these two functions for a given problem will produce the
fastest possible code.


## Examples

An example script that simulates decaying two-dimensional turbulence reproducing
the results of the paper by

  > McWilliams, J. C. (1984). The emergence of isolated coherent vortices in turbulent flow. _J. Fluid Mech._, **146**, 21-43.

is found in `examples/twodturb/McWilliams.jl`.



## Future work

The code is in the chaos stage of development. A main goal for the future
is to permit the use of shared memory parallelism in the compute-intensive
routines (shared-memory parallelism provided already by FFTW/MKLFFT, but
is not yet native to Julia for things like element-wise matrix multiplication,
addition, and assignment). This feature may possibly be enabled by
Intel Lab's
[ParallelAccelerator](https://github.com/IntelLabs/ParallelAccelerator.jl)
package.

## Authors

FourierFlows is currently being developed by [Gregory L. Wagner](https://glwagner.github.io) and [Navid C. Constantinou](http://www.navidconstantinou.com).

## DocStrings

```@contents
Pages = [
    "man/docstrings.md",
    ]
Depth = 2
```

## Index

```@index
Pages = [
    "man/docstrings.md",
    ]
```