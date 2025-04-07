using MPI
@info "Initialising MPI job"
MPI.Init()


using BenchmarkTools
using LinearAlgebra
using Logging
using Parameters
using WAVI

include("grid_dd.jl");

function run_mpi()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    comm_size = MPI.Comm_size(comm)

    if rank == 0
        @info "Starting job on rank $(rank)/$(comm_size)"
        @warn "This implementation does not work in any way, please do not use"
        p = GridParams(100, 100; mx=2, my=2)
        mpi_spec = MPISpec(
            ngridsx = p.mx,
            ngridsy = p.my,
            niterations = 10
        )
        model = create_model(p, mpi_spec)
        run_grid_ops(model)
    end
end

run_mpi()
