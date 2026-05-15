# WAVI Model Specifications

## Overview

A late addition to WAVI is the idea of specifications, which guide the setup of the library. 

Semantically, we can think of the components of WAVI this way: 

* `Model` brings together all the elements of compute and data into a singular operable entity
* `Simulation` determines how the model will be used to produce results over time
* `Specs` determines how the model and simulation work together to solve the ice sheet system

Specifications heavily leverage multiple dispatch to re-implement portions of WAVI to introduce new structural architectures, computational approaches and data operations.

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

## Domain Decomposition

### ThreadedSpec

The threaded specification utilises an iterative Schwarz approach to domain decomposition, producing individual domains prior to execution of the velocity solve processing. 

For this specification only two methods are overridden:

`update_preconditioner!` creates a decomposed grid across individual threads, allowing these threads to `precondition!` individually within that subdomain.

`precondition!` then handles exchange of velocities takes place during each iteration of the velocity solve. The state is updated following a transfer of velocities from global to the "local" grid, upon which the local `update_state!` is used to update the local model prior to transferring back to the global grid using a [partition of unity](https://en.wikipedia.org/wiki/Partition_of_unity).

The memory space is still limited to one process, with the potential to leverage threading optimisations within a single physical processor. _Therefore this might increase parallelisation of the model computations, __but it will remain memory bound__._

### MPISpec

The MPI specification is used to decompose the global domain of the `Model`. Unlike the threaded specification however, this decomposition takes place at time of creation. The grid provided is split and the root rank takes the upper-left node in the topology, [which is clearly explained here](https://hpc-tutorials.llnl.gov/mpi/virtual_topologies/).

The list of methods overridden is not explained in detail here, there are numerous areas of the model that require alteration:

* The underlying `Model` constructor is implemented such that each node contains it's own localised portion of the global domain
* Velocity solves incorporate a halo exchange after the computation is carried out
* `Simulation` operations are updated to collect data as required for outputs using a deferred evaluation mechanism, required to ensure outputs are taken from the global grid and rendered from the root process only

The topology of an individual ranks grid (in this case rank zero) is also relevant, with velocities at the neighboured-edges fixed and then transferred, in much the same manner as with ThreadedSpec. The following diagram illustrates the situation:

```@raw html
<center><img src="../assets/mpispec_exchange.png" alt="" title="" height="700" /></center>
```

#### Distributed TODO (MPISpec)

There are elements of MPISpec still under development at time of writing: 

* [Damping and partitions of unity are not yet applied](https://github.com/WAVI-ice-sheet-model/WAVI.jl/issues/108)
* [Outputs only can be taken from `model.global_fields`](https://github.com/WAVI-ice-sheet-model/WAVI.jl/issues/111) - therefore melt rates, parameters and other fields are not accessible as outputs
* Non-2D fields are not accessible as distributed outputs

## Specifications API

Please refer to the [API documentation for the structure / constructor definitions](API/specifications.md).