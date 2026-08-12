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
* **Grid Balancing**: If the global grid size is not perfectly divisible by the number of subdomains or MPI processes, the remainder cells are balanced across the first few subdomains (for both `ThreadedSpec` and `MPISpec`). This ensures the full global domain is covered without remainders at the boundaries, removing the previous requirement of the user needing to specify a divisible thread/core count.

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

Configure MPI for your Julia project once per machine (MPIPreferences, `mpiexecjl`, verification). See [MPI setup](./mpi_setup.md) and the [MPI.jl configuration guide](https://juliaparallel.org/MPI.jl/stable/configuration/).

To run a driver with an `MPISpec`:

```bash
mpiexecjl -n <num> --project=<path-to-project> julia <path-to-driver.jl> <driver args...>
```

Example (MISMIP+ driver, four ranks):

```bash
mpiexecjl -n 4 --project=../.. julia example_drivers/MISMIP_PLUS/MISMIP_PLUS.jl
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

* The underlying `Model` constructor is implemented such that each node contains it's own localised portion of the global domain.

!!! tip "Process Layout Geometry"
    The geometric layout of the MPI process grid (`px` × `py`) heavily impacts performance and solver convergence. For domains with high aspect ratios (like MISMIP+ which is long and narrow), **favour 1D process layouts** (e.g., `4x1` instead of `2x2`).

    A `2x2` layout on a narrow grid creates subdomains with very few "core" cells compared to the halo size, which drastically increases the relative overhead of Schwarz Partition-of-Unity (`pou`) communication and slows down convergence.

* Velocity solves use overlapping **Schwarz iterations** (multiple local subdomain solves per step). The subdomain boundary treatment is controlled by the `pou` parameter:
    * **`pou=true` (default)**: Uses a **Partition of Unity** approach. Subdomain velocities are smoothly blended together in overlap zones, ensuring physical continuity. This requires a `halo` size of at least `2`.
    * **`pou=false`**: Uses **Restricted Additive Schwarz (RAS)**. Velocities are exchanged directly between neighboring subdomains without blending. This is computationally simpler and allows a smaller `halo` size of `1`.
* **Thickness Halos**: Ice thickness `h` is synchronised across subdomains at the end of each timestep.
* **Global Residual Check**: Picard convergence is verified globally across all ranks. Each process computes the squared norm of its local residual on "core" unknowns only (excluding halo overlaps) to prevent double-counting. These are then aggregated via `MPI.Allreduce` to determine the global relative residual.
* `Simulation` operations are updated to collect data as required for outputs using a deferred evaluation mechanism, required to ensure outputs are taken from the global grid and rendered from the root process only.

The topology of an individual ranks grid (in this case rank zero) is also relevant, with velocities at the neighboured-edges fixed and then transferred. The following diagram illustrates the situation:

```@raw html
<center><img src="../assets/mpispec_exchange.png" alt="" title="" height="700" /></center>
```

#### Distributed TODO (MPISpec)

There are elements of MPISpec still under development at time of writing: 

* [Outputs only can be taken from `model.global_fields`](https://github.com/WAVI-ice-sheet-model/WAVI.jl/issues/111) - therefore melt rates, parameters and other fields are not accessible as outputs (though `model.mpi_rank` is now available).
* Non-2D fields (e.g., 3D temperature) are not yet accessible as distributed outputs.
* Checkpointing writes a portable full-domain `BasicSpec` snapshot (restart fields); see [Checkpoints](./API/timestepping_params.md#mpispec-checkpoints).

## Specifications API

Please refer to the [API documentation for the structure / constructor definitions](API/specifications.md).