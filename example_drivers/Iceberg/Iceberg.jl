# Derived from test/verification_tests/iceberg_test.jl
using MPI
using WAVI

function Iceberg_Grid(;
        nx::Int = 6,
        ny::Int = 4,
    )::Grid
    nσ = 4
    x0 = -50000.0
    y0 = -50000.0
    dx = 5000.0
    dy = 5000.0

    # Inhomogenous Dirichlet boundary conditions will be applied on x- and y- center lines
    u_isfixed = falses(nx+1, ny)
    u_isfixed[div(nx,2)+1, :] .= true
    v_isfixed = falses(nx, ny+1)
    v_isfixed[:, div(ny,2)+1] .= true

    grid = Grid(nx = nx, 
                ny = ny, 
                nσ = nσ, 
                x0 = x0, 
                y0 = y0, 
                dx = dx, 
                dy = dy,
                u_isfixed = u_isfixed,
                v_isfixed = v_isfixed)

    return grid
end

function Iceberg(;
        folder = "outputs",
        grid = Iceberg_Grid(),
        spec = BasicSpec(),
        end_time = 1000.
    )
    
    # Physical parameters & Initial Conditions
    nx, ny = grid.nx, grid.ny
    
    # # Model domain mask on h-grid
    starting_thickness = zeros(nx, ny) .+ 200.0
    initial_conditions = InitialConditions(initial_thickness = starting_thickness)

    accumulation_rate = 0.3
    default_temperature = 265.700709
    params = Params(accumulation_rate = accumulation_rate,
                    default_temperature = default_temperature)

    bed_elevation = -500.0 .* ones(nx, ny)

    maxiter_picard = 1 # No need for Picard iteration if running to steady state
    solver_params = SolverParams(maxiter_picard = maxiter_picard)

    # Make the model
    model = Model(grid = grid,
                  bed_elevation = bed_elevation,
                  params = params,
                  solver_params = solver_params,
                  initial_conditions = initial_conditions,
                  spec = spec)

    # Timestepping
    chkpt_freq = 20.0
    timestepping_params = TimesteppingParams(dt = 0.1, 
                                             end_time = end_time,
                                             chkpt_freq = chkpt_freq)

    # Output parameters (following MISMIP_PLUS example)
    outputs = (
        h = "global_fields.gh.h",
        u = "global_fields.gh.u",
        v = "global_fields.gh.v",
    )
    # Use a reasonable output frequency or defaults
    output_freq = 20. 
    output_params = OutputParams(outputs,
                                 output_path = folder,
                                 output_freq = output_freq,
                                 output_format = "jld2",
                                 zip_format = "nc",
                                 output_start = true,
                                )

    simulation = Simulation(model = model, 
                            timestepping_params = timestepping_params,
                            output_params = output_params)
    
    run_simulation!(simulation)
    return simulation
end

if abspath(PROGRAM_FILE) == @__FILE__
    MPI.Init()
    if MPI.Comm_size(MPI.COMM_WORLD) > 1
        # MPI Mode
        grid = Iceberg_Grid()
        mpi_spec = MPISpec(MPI.Comm_size(MPI.COMM_WORLD), 1, 2, grid; pou=true, niterations=2)

        Iceberg(
            folder = "outputs/iceberg_mpi",
            grid = grid,
            spec = mpi_spec
        )
    elseif Threads.nthreads() > 1
        """
        Threaded Mode
        Run with:
        julia --project -t 2 example_drivers/Iceberg/Iceberg.jl
        """
        grid = Iceberg_Grid()
        threaded_spec = ThreadedSpec(ngridsx=Threads.nthreads(), ngridsy=1, overlap=2, niterations=2)
        Iceberg(
            folder = "outputs/iceberg_thread",
            grid = grid,
            spec = threaded_spec,
        )
    else
        # Serial Mode
        Iceberg(folder = "outputs/iceberg_serial")
    end
end
