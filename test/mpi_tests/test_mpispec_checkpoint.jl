# Unit test: MPISpec checkpoint write + pickup (run under MPI with 2 ranks).
using Test
using MPI
using WAVI

@testset "MPISpec checkpoint write/pickup" begin
    MPI.Initialized() || MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)
    nprocs == 2 || error("MPISpec checkpoint test expects 2 ranks, got $nprocs")

    nx, ny, nσ = 12, 8, 3
    grid = Grid(nx = nx, ny = ny, nσ = nσ, dx = 1.0, dy = 1.0)
    bed = -500.0 .* ones(nx, ny)
    ics = InitialConditions(
        initial_thickness = 100.0 .* ones(nx, ny),
        initial_temperature = fill(263.15, nx, ny, nσ),
        initial_viscosity = fill(1.0e15, nx, ny, nσ),
    )
    params = Params(accumulation_rate = 0.1)
    solver_params = SolverParams(maxiter_picard = 1)

    # Shared path across ranks (getpid differs per process).
    chkpt_dir = joinpath(tempdir(), "wavi_mpi_chkpt_test_12x8")
    rank == 0 && (isdir(chkpt_dir) && rm(chkpt_dir; force = true, recursive = true))
    MPI.Barrier(comm)
    rank == 0 && mkpath(chkpt_dir)
    MPI.Barrier(comm)

    function build_model()
        spec = MPISpec(nprocs, 1, 2, grid; pou = true, niterations = 2)
        return Model(
            grid = grid,
            bed_elevation = bed,
            params = params,
            solver_params = solver_params,
            initial_conditions = ics,
            spec = spec,
            verbose = false,
        )
    end

    function gather_hθη(model)
        return (
            copy(WAVI.Specs.collect_mpi_field!(model, [:global_fields, :gh, :h])),
            copy(WAVI.Specs.collect_mpi_field!(model, [:global_fields, :g3d, :θ])),
            copy(WAVI.Specs.collect_mpi_field!(model, [:global_fields, :g3d, :η])),
        )
    end

    # Run to a checkpoint boundary so end state matches the file used for pickup.
    model1 = build_model()
    tp1 = TimesteppingParams(
        dt = 0.1,
        end_time = 0.2,
        chkpt_freq = 0.2,
        chkpt_path = chkpt_dir,
    )
    sim1 = Simulation(
        model = model1,
        timestepping_params = tp1,
        output_params = OutputParams(output_path = chkpt_dir),
    )
    run_simulation!(sim1)

    h_after_first, θ_after_first, η_after_first = gather_hθη(sim1.model)
    n_iter_chkpt = sim1.timestepping_params.n_iter_chkpt
    @test n_iter_chkpt > 0

    if rank == 0
        fname = joinpath(chkpt_dir, WAVI.Outputs.checkpoint_filename(n_iter_chkpt))
        @test isfile(fname)
    end
    MPI.Barrier(comm)

    # Continuous control to the same final time as the restarted run.
    model_ctrl = build_model()
    tp_ctrl = TimesteppingParams(dt = 0.1, end_time = 0.5, chkpt_freq = Inf)
    sim_ctrl = Simulation(
        model = model_ctrl,
        timestepping_params = tp_ctrl,
        output_params = OutputParams(),
    )
    run_simulation!(sim_ctrl)
    h_ctrl, θ_ctrl, η_ctrl = gather_hθη(sim_ctrl.model)

    # Pickup from checkpoint and continue to the same end time.
    model2 = build_model()
    tp2 = TimesteppingParams(
        niter0 = n_iter_chkpt,
        dt = 0.1,
        end_time = 0.5,
        chkpt_freq = Inf,
        chkpt_path = chkpt_dir,
    )
    sim2 = Simulation(
        model = model2,
        timestepping_params = tp2,
        output_params = OutputParams(output_path = chkpt_dir),
    )
    @test sim2.clock.n_iter == n_iter_chkpt
    @test sim2.model.spec isa MPISpec

    # Immediately after pickup: restored fields match the checkpointed first-run state.
    h_pickup, θ_pickup, η_pickup = gather_hθη(sim2.model)
    if rank == 0
        @test all(isfinite.(h_pickup))
        @test all(isfinite.(θ_pickup))
        @test all(isfinite.(η_pickup))
        @test isapprox(h_pickup, h_after_first; rtol = 1e-5, atol = 1e-6)
        @test isapprox(θ_pickup, θ_after_first; rtol = 1e-5, atol = 1e-6)
        @test isapprox(η_pickup, η_after_first; rtol = 1e-5, atol = 1e-6)
    end
    MPI.Barrier(comm)

    run_simulation!(sim2)
    h_restart, θ_restart, η_restart = gather_hθη(sim2.model)

    if rank == 0
        @test all(isfinite.(h_after_first))
        @test all(isfinite.(θ_after_first))
        @test all(isfinite.(η_after_first))
        # Restarted run progressed beyond the pre-restart state.
        @test !isapprox(h_restart, h_after_first; rtol = 1e-5, atol = 1e-6)
        @test isapprox(h_restart, h_ctrl; rtol = 1e-5, atol = 1e-6)
        @test isapprox(θ_restart, θ_ctrl; rtol = 1e-5, atol = 1e-6)
        @test isapprox(η_restart, η_ctrl; rtol = 1e-5, atol = 1e-6)
        rm(chkpt_dir; force = true, recursive = true)
    end
    MPI.Barrier(comm)
end

MPI.Finalize()
