using MPI
using WAVI
using WAVI.Specs

export @root

"""
    Get directional halos for a given MPISpec

    This function returns the halos for a given MPISpec in the form of a tuple (top, right, bottom, left).
    If a neighbour is not present, the halo is set to 0.

    Args:
        spec: The MPISpec to get the halos for

    Returns:
        A tuple of the form (top, right, bottom, left)
"""
function get_halos(spec::MPISpec)::Tuple{Int, Int, Int, Int}
    # Note: MPI_Cart_shift returns `MPI_PROC_NULL` (represented as a negative integer)
    #       if the neighbour is not present
    return spec.top > -1 ? spec.halo : 0, 
           spec.right > -1 ? spec.halo : 0, 
           spec.bottom > -1 ? spec.halo : 0, 
           spec.left > -1 ? spec.halo : 0
end

"""
    Get the bounds of a given MPISpec

    The bounds are calculated based on the global grid size and the number of processes in each direction.

    Args:
        spec: The MPISpec to get the bounds for

    Returns:
        A tuple of the form (x_start, x_end, y_start, y_end)
"""
function get_bounds(spec::MPISpec)::Tuple{Int, Int, Int, Int}
    return max(spec.coords[1] * div(spec.global_grid.nx, spec.px) + 1 - get_halos(spec)[4], 1),
           min(spec.coords[1] * div(spec.global_grid.nx, spec.px) + div(spec.global_grid.nx, spec.px) + get_halos(spec)[2], spec.global_grid.nx),
           max(spec.coords[2] * div(spec.global_grid.ny, spec.py) + 1 - get_halos(spec)[1], 1),
           min(spec.coords[2] * div(spec.global_grid.ny, spec.py) + div(spec.global_grid.ny, spec.py) + get_halos(spec)[3], spec.global_grid.ny)
end


"""
    Get the size of a given MPISpec

    The size is calculated based on the global grid size and the number of processes in each direction.

    Args:
        spec: The MPISpec to get the size for

    Returns:
        A tuple of the form (nx, ny)
"""
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