using MPI
using WAVI 

include("MISMIP_PLUS.jl")
include("plotting.jl")
include("utils.jl")

if abspath(PROGRAM_FILE) == @__FILE__
    grid = MISMIP_PLUS_GRID()
    nprocs = 2          # TODO: Retrieve from mpiexecjl invocation arguments
    mpi_spec = MPISpec(nprocs, 1, 2, grid)

    benchmark_main("mpi.$(nprocs)", MISMIP_PLUS, Dict{Symbol, Any}(
        :grid => grid,
        :spec => mpi_spec,
    ), ["h", "u", "v"], mpi_spec.rank)
    MPI.Finalize()
end