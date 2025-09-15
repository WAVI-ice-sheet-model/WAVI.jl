# WAVI Model Specifications

## Overview

A late addition to WAVI is the idea of "Specs", which guide the setup of the model itself. These can be seen to operate underneath `Simulation`, which brings together `Model` and `Clock` to apply a lifecycle for WAVI as an full model. 

Semantically, we can think of the components this way: 

* `Model` brings together all the elements of compute and data into a singular operable entity
* `Spec` determines how the model actually operates at the computational level
* `Simulation` determines how the model will be used to produce results over time

## Types

### Model Specifications

There are three types of specifications currently available in WAVI: 

* `BasicSpec` - this is the original single-core, single-thread implementation. You will run into limits with bigger domains.
* `ThreadedSpec` - this is an implementation of schwarz decomposition for solving subdomains within individual threads. It won't get you around the memory limits, but it will speed up compute.
* `MPISpec` - this is the __currently experimental__ fully-distributed specification for WAVI. The model domain is split apart and computed individually.

### Implementation Considerations

Depending on specifications, some considerations need to be taken into account:

* Examples in the docs simply default to `BasicSpec`, so the model will not be distributed by default
* How you use `Model` interface changes in a distributed setting, to a fully functional interface

## Setting up specifications

This is relatively simply achieved in the call to Model. As part of the model creation process, one needs to supply a [`Grid`](@ref) to a [`Model`](@ref) constructor (alongside the bed). If no spec is provided a `BasicSpec` is used.

For any other specification, you might use:

```julia
    spec = ThreadedSpec(ngridsx = 16, 
                        ngridsy = 2,
                        overlap = 1,
                        niterations = 1)

    # or

    # px, py, halo, grid
    spec = MPISpec(16, 2, 1, grid)
```

As can be seen, the grid is necessarily provided to the MPI specification prior to model construction, whereas the operation of the threaded specification does not require the same.

#### MPI specific execution

MPI.jl has a dependency on MPIPreferences in order to set up the necessary links to your execution framework. This is outside the scope here to pick up, but you [should follow these instructions](https://juliaparallel.org/MPI.jl/stable/usage/#Julia-wrapper-for-mpiexec) in order to access `mpiexecjl`.

To then run with an `MPISpec`, use the following command: 

```julia
mpiexecjl --project -np <num> julia <path to driver> <driver args...>`
```

## Specifications API

```@docs; canonical=false
BasicSpec
ThreadedSpec
MPISpec
```