using MPI
using WAVI
using WAVI.Specs

export @root

function get_halos(spec::MPISpec)::Tuple{Int, Int, Int, Int}
    return spec.top > -1 ? spec.halo : 0, 
           spec.right > -1 ? spec.halo : 0, 
           spec.bottom > -1 ? spec.halo : 0, 
           spec.left > -1 ? spec.halo : 0
end

function get_bounds(spec::MPISpec)::Tuple{Int, Int, Int, Int}
    return max(spec.coords[1] * div(spec.global_grid.nx, spec.px) + 1 - get_halos(spec)[4], 1),
           min(spec.coords[1] * div(spec.global_grid.nx, spec.px) + div(spec.global_grid.nx, spec.px) + get_halos(spec)[2], spec.global_grid.nx),
           max(spec.coords[2] * div(spec.global_grid.ny, spec.py) + 1 - get_halos(spec)[1], 1),
           min(spec.coords[2] * div(spec.global_grid.ny, spec.py) + div(spec.global_grid.ny, spec.py) + get_halos(spec)[3], spec.global_grid.ny)
end

function get_size(spec::MPISpec)::Tuple{Int, Int}
    return div(spec.global_grid.nx, spec.px) + get_halos(spec)[4] + get_halos(spec)[2], 
           div(spec.global_grid.ny, spec.py) + get_halos(spec)[1] + get_halos(spec)[3]
end

# Sourced from Oceananigans.DistributedComputations
# - available under the MIT License: https://github.com/CliMA/Oceananigans.jl/blob/main/LICENSE
# - original copyright of the below code (except adaptations for use) 2018 Climate Modeling Alliance

mpi_initialized()     = MPI.Initialized()
mpi_rank(comm)        = MPI.Comm_rank(comm)
mpi_size(comm)        = MPI.Comm_size(comm)
global_barrier(comm)  = MPI.Barrier(comm)
global_communicator() = MPI.COMM_WORLD

"""
    @root communicator exp...

Perform `exp` only on rank 0 in communicator, otherwise known as the "root" rank.
Other ranks will wait for the root rank to finish before continuing.
If `communicator` is not provided, `MPI.COMM_WORLD` is used.
"""
macro root(communicator, exp)
    command = quote
        if WAVI.Specs.mpi_initialized()
            rank = WAVI.Specs.mpi_rank($communicator)
            if rank == 0
                $exp
            end
            WAVI.Specs.global_barrier($communicator)
        else
            $exp
        end
    end
    return esc(command)
end

macro root(exp)
    command = quote
        @root WAVI.Specs.global_communicator() $exp
    end
    return esc(command)
end

## END Sourced from Oceananigans.DistributedComputations