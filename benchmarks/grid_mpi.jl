using MPI
@info "Initialising MPI job"
MPI.Init()


using BenchmarkTools
using LinearAlgebra
using Logging
using Parameters
using WAVI

include("grid_dd.jl");


comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm) + 1
comm_size = MPI.Comm_size(comm)

if rank == 1
    @info "Initialising full domain in rank $(rank)/$(comm_size)"

else
    @info "Initialising sub domain model in rank $(rank)/$(comm_size)"

end
