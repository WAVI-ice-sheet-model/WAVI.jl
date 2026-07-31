using PrecompileTools
using MPI
using Base.CoreLogging: with_logger, NullLogger

@setup_workload begin
    # Create a tiny 10x10 model.
    # Just big enough to force Julia to compile the math, but
    # small enough that it should compiles quickly without
    # slowing down installation.
    grid = Grid(
        nx = 10, ny = 10, nσ = 2,
        x0 = 0.0, y0 = -40000.0,
        dx = 8000.0, dy = 8000.0,
        h_mask = trues(10, 10),
        u_iszero = falses(11, 10),
        v_iszero = falses(10, 11),
    )
    bed = zeros(10, 10)
    params = Params()
    solver_params = SolverParams(maxiter_picard = 1)

    @compile_workload begin
        # Run a standard model (BasicSpec). This forces Julia to compile a large amount
        # of the core math and physics that every run needs.
        model_basic = Model(grid, bed, BasicSpec(); params = params, solver_params = solver_params)
        update_velocities!(model_basic)

        # Next, run the threaded model.
        try
            spec_threaded = ThreadedSpec(; ngridsx = 2, ngridsy = 2, overlap = 1, niterations = 1)
            model_threaded = Model(grid, bed, spec_threaded; params = params, solver_params = solver_params)
            update_velocities!(model_threaded)
        catch
            # Do not crash the installation if the grid is too small for threaded decomposition.
        end

        # Finally, try to run a tiny MPI model on a single processor to
        # compile the remaining MPI-specific scatter/gather code. However, because
        # Pkg.precompile often crashes when it tries to start MPI on some servers, this is
        # wrapped in a safe try...catch. If it fails, it just silently gives up, but hopefully
        # majority of the speedup is still gained from the blocks above.
        if get(ENV, "WAVI_PRECOMPILE_MPI", "1") != "0"
            try
                with_logger(NullLogger()) do
                    MPI.Initialized() || MPI.Init()
                    if MPI.Comm_size(MPI.COMM_WORLD) == 1
                        # Force the halo exchange to run even on a 1x1 grid just to compile
                        # the MPI methods.
                        spec_mpi = MPISpec(1, 1, 1, grid; pou = false, niterations = 1)
                        model_mpi = Model(
                            grid, bed, spec_mpi;
                            params = params,
                            solver_params = solver_params,
                            verbose = false,
                        )
                        update_velocities!(model_mpi)
                        Specs.halo_exchange!(model_mpi; fields = [:h, :u, :v])
                        Specs.collect_mpi_field!(model_mpi, [:global_fields, :gh, :h])
                    end
                end
            catch
                # MPI crashed or is unavailable. Silently ignore the error so the installation
                # completes safely.
            end
        end
    end
end
