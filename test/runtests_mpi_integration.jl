using MPI
using Test

include("run_mpi_script.jl")

@testset "MPISpec MPI integration tests" begin
    p, cmd = run_mpi_script("test_mpispec_iceberg_integration.jl"; nprocs = 2)
    @test success(p)
    if !success(p)
        @error "MPI integration tests failed" cmd exitcode = p.exitcode
    end
end
