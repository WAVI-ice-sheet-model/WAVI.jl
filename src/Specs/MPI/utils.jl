using MPI
using WAVI
using WAVI.Specs

export @root

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