using MPI
using Test

include("run_mpi_script.jl")

@testset "MPISpec MPI unit tests" begin
    p, cmd = run_mpi_script("test_mpispec_halo_collect.jl"; nprocs = 2)
    @test success(p)
    if !success(p)
        @error "MPI unit tests failed" cmd exitcode = p.exitcode
    end

    p, cmd = run_mpi_script("test_mpispec_checkpoint.jl"; nprocs = 2)
    @test success(p)
    if !success(p)
        @error "MPI checkpoint tests failed" cmd exitcode = p.exitcode
    end
end
