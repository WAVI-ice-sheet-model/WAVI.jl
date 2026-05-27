# Diagnostic (not run by CI): measure max relative error MPI vs BasicSpec for Iceberg setups.
# Run from repo root:
#   mpiexecjl -n 2 --project=test julia test/mpi_tests/measure_iceberg_mpi_error.jl
# Use the output to tune _MPI_ICEBERG_RTOL in test_mpispec_iceberg_integration.jl.
using MPI
using WAVI
using Statistics

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

function build_simulation(spec; end_time = 2.0, maxiter_picard = 1)
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
    tp = TimesteppingParams(dt = 0.1, end_time = end_time, chkpt_freq = Inf)
    sim = Simulation(model = model, timestepping_params = tp, output_params = OutputParams())
    run_simulation!(sim)
    return sim
end

"""
    field_errors(a, b)

Compare two fields and return max absolute error, max relative error, and mean relative error.
"""
function field_errors(a, b)
    diff = abs.(a .- b)
    denom = max.(abs.(b), 1e-12)
    rel = diff ./ denom
    return (
        max_abs = maximum(diff),
        max_rel = maximum(rel),
        mean_rel = mean(rel),
    )
end

"""
    main()

Print MPI vs BasicSpec field errors for several solver settings (rank 0 only).
"""
function main()
    MPI.Initialized() || MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)

    configs = [
        (; label = "baseline", pou = false, niterations = 2, maxiter_picard = 1),
        (; label = "niter5", pou = false, niterations = 5, maxiter_picard = 1),
        (; label = "picard5", pou = false, niterations = 5, maxiter_picard = 5),
        (; label = "picard10", pou = false, niterations = 5, maxiter_picard = 10),
        (; label = "pou_picard5", pou = true, niterations = 5, maxiter_picard = 5),
    ]

    grid = iceberg_grid()

    for cfg in configs
        if rank == 0
            basic = build_simulation(BasicSpec(); maxiter_picard = cfg.maxiter_picard)
            h_basic = basic.model.fields.gh.h
            u_basic = basic.model.fields.gu.u
            v_basic = basic.model.fields.gv.v
        end
        MPI.Barrier(comm)
        mpi_spec = MPISpec(nprocs, 1, 2, grid; pou = cfg.pou, niterations = cfg.niterations)
        mpi_sim = build_simulation(mpi_spec; maxiter_picard = cfg.maxiter_picard)
        h_mpi = WAVI.Specs.collect_mpi_field!(mpi_sim.model, [:global_fields, :gh, :h])
        u_mpi = WAVI.Specs.collect_mpi_field!(mpi_sim.model, [:global_fields, :gu, :u])
        v_mpi = WAVI.Specs.collect_mpi_field!(mpi_sim.model, [:global_fields, :gv, :v])
        if rank == 0
            println("\n=== $(cfg.label) (pou=$(cfg.pou), niter=$(cfg.niterations), picard=$(cfg.maxiter_picard)) ===")
            for (name, a, b) in [(:h, h_mpi, h_basic), (:u, u_mpi, u_basic), (:v, v_mpi, v_basic)]
                e = field_errors(a, b)
                println("$name: max_rel=$(e.max_rel) max_abs=$(e.max_abs) mean_rel=$(e.mean_rel)")
            end
        end
        MPI.Barrier(comm)
    end
end

main()
MPI.Finalize()
