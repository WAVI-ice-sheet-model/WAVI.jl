# Integration test: short Iceberg run on MPISpec should match BasicSpec (run under MPI with 2 ranks).
using Test
using MPI
using WAVI

# Tolerances below are set from measured max relative error (MPI vs BasicSpec) for this
# setup with pou=true, niterations=5, maxiter_picard=5, end_time=2.0, 2 ranks:
#   h max_rel ≈ 7.8e-5, u ≈ 6.3e-4, v ≈ 8.1e-4 (see measure_iceberg_mpi_error.jl).
const _MPI_ICEBERG_RTOL = Dict(:h => 2e-4, :u => 2e-3, :v => 2e-3)

# Short Iceberg run: gathered MPISpec fields on rank 0 should match a serial BasicSpec run.
@testset "MPISpec Iceberg integration comparison" begin
    MPI.Initialized() || MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)

    @test nprocs == 2

    """
        iceberg_grid(; nx, ny)

    Build the small Iceberg test grid with fixed velocities on the centre lines.
    """
    function iceberg_grid(; nx = 12, ny = 8)
        u_isfixed = falses(nx + 1, ny)
        u_isfixed[div(nx, 2) + 1, :] .= true
        v_isfixed = falses(nx, ny + 1)
        v_isfixed[:, div(ny, 2) + 1] .= true
        return Grid(
            nx = nx,
            ny = ny,
            nσ = 4,
            x0 = -50000.0,
            y0 = -50000.0,
            dx = 5000.0,
            dy = 5000.0,
            u_isfixed = u_isfixed,
            v_isfixed = v_isfixed,
        )
    end

    """
        build_simulation(spec; end_time, maxiter_picard)

    Run a short Iceberg simulation with the given model spec and return the simulation.
    """
    function build_simulation(spec; end_time = 2.0, maxiter_picard = 5, chkpt_freq = Inf, niter0 = 0)
        grid = iceberg_grid()
        nx, ny = grid.nx, grid.ny
        model = Model(
            grid = grid,
            bed_elevation = -500.0 .* ones(nx, ny),
            params = Params(accumulation_rate = 0.3, default_temperature = 265.700709),
            solver_params = SolverParams(maxiter_picard = maxiter_picard),
            initial_conditions = InitialConditions(initial_thickness = 200.0 .* ones(nx, ny)),
            spec = spec,
        )
        tp = TimesteppingParams(dt = 0.1, end_time = end_time, chkpt_freq = chkpt_freq, niter0 = niter0)
        sim = Simulation(model = model, timestepping_params = tp, output_params = OutputParams())
        run_simulation!(sim)
        return sim
    end

    mktempdir() do temporary_directory
      cd(temporary_directory) do
        println("Running MPI integration and checkpoint test in directory:",temporary_directory)

        grid = iceberg_grid()
        # Match MPISpec defaults used in Iceberg MPI drivers: PoU + default Schwarz iterations.
        mpi_spec = MPISpec(nprocs, 1, 2, grid; pou = true, niterations = 5)
        mpi_sim = build_simulation(mpi_spec; chkpt_freq = 1.0)
        mpi_sim_pickup = build_simulation(mpi_spec; chkpt_freq = 1.0, niter0=10)

        h_mpi = WAVI.Specs.collect_mpi_field!(mpi_sim.model, [:global_fields, :gh, :h])
        u_mpi = WAVI.Specs.collect_mpi_field!(mpi_sim.model, [:global_fields, :gu, :u])
        v_mpi = WAVI.Specs.collect_mpi_field!(mpi_sim.model, [:global_fields, :gv, :v])

        h_mpi_pickup = WAVI.Specs.collect_mpi_field!(mpi_sim_pickup.model, [:global_fields, :gh, :h])
        u_mpi_pickup = WAVI.Specs.collect_mpi_field!(mpi_sim_pickup.model, [:global_fields, :gu, :u])
        v_mpi_pickup = WAVI.Specs.collect_mpi_field!(mpi_sim_pickup.model, [:global_fields, :gv, :v])

        if rank == 0
            basic_sim = build_simulation(BasicSpec())
            h_basic = basic_sim.model.fields.gh.h
            u_basic = basic_sim.model.fields.gu.u
            v_basic = basic_sim.model.fields.gv.v

            @test all(isfinite.(h_mpi))
            @test all(isfinite.(u_mpi))
            @test all(isfinite.(v_mpi))

            @test isapprox(h_mpi, h_basic; rtol = _MPI_ICEBERG_RTOL[:h], atol = 1e-6)
            @test isapprox(u_mpi, u_basic; rtol = _MPI_ICEBERG_RTOL[:u], atol = 1e-6)
            @test isapprox(v_mpi, v_basic; rtol = _MPI_ICEBERG_RTOL[:v], atol = 1e-6)

            @test all(isfinite.(h_mpi_pickup))
            @test all(isfinite.(u_mpi_pickup))
            @test all(isfinite.(v_mpi_pickup))

            @test h_mpi_pickup == h_mpi
            @test u_mpi_pickup == u_mpi
            @test v_mpi_pickup == v_mpi

        end
      end
    end
end

MPI.Finalize()
