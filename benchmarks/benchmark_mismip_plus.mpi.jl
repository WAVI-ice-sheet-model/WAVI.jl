using MPI
using WAVI 

include("MISMIP_PLUS.jl")
include("plotting.jl")
include("utils.jl")

if abspath(PROGRAM_FILE) == @__FILE__
    grid = MISMIP_PLUS_GRID()
    mpi_spec = MPISpec(2, 1, 2, grid)

    benchmark_main("mpi", MISMIP_PLUS, Dict{Symbol, Any}(
        :grid => grid,
        :spec => mpi_spec,
    ), ["h", "u", "v"], mpi_spec.rank)
end