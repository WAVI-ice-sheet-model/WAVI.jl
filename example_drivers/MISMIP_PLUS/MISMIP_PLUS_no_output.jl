using MPI
using WAVI 

function MISMIP_PLUS()
    #Grid and boundary conditions
    nx = 160
    ny = 20
    nσ = 4
    x0 = 0.0
    y0 = -40000.0
    dx = 4000.0
    dy = 4000.0
    h_mask=trues(nx,ny)
    u_iszero = falses(nx+1,ny); u_iszero[1,:].=true
    v_iszero=falses(nx,ny+1); v_iszero[:,1].=true; v_iszero[:,end].=true
    grid = Grid(nx = nx, 
                ny = ny,   
                nσ = nσ, 
                x0 = x0, 
                y0 = y0, 
                dx = dx, 
                dy = dy,
                h_mask = h_mask, 
                u_iszero = u_iszero, 
                v_iszero = v_iszero)

    #Bed 
    bed = WAVI.mismip_plus_bed #function definition

    #solver parameters
    maxiter_picard = 1
    super_implicitness = 1.0
    solver_params = SolverParams(maxiter_picard = maxiter_picard, super_implicitness=super_implicitness)

    #spec = BasicSpec()
    spec = ThreadedSpec(ngridsx = 16,ngridsy=2,overlap=1,niterations=1)
    #spec = MPISpec(16, 2, 1, grid)

    #Physical parameters
    default_thickness = 100.0 #set the initial condition this way
    accumulation_rate = 0.3
    params = Params(default_thickness = default_thickness, 
                    accumulation_rate = accumulation_rate)

    #make the model
    model = Model(grid,
                  bed,
                  spec; 
                  params = params, 
                  solver_params = solver_params)
 
    #timestepping parameters
    niter0 = 0
    dt = 0.1
    end_time = 10000.0
    timestepping_params = TimesteppingParams(niter0 = niter0, 
                                             dt = dt, 
                                             end_time = end_time)
    
    simulation = Simulation(model = model, timestepping_params = timestepping_params)
            
    #perform the simulation
    run_simulation!(simulation)
    return simulation
end

if abspath(PROGRAM_FILE) == @__FILE__
    # A little bootstapping way of running the MISMIP+ experiment
    MPI.Init()
    if MPI.Comm_size(MPI.COMM_WORLD) > 1 
        grid = MISMIP_PLUS_GRID()
        mpi_spec = MPISpec(MPI.Comm_size(MPI.COMM_WORLD), 1, 2, grid)

        sim = @time MISMIP_PLUS(
            grid = grid,
            spec = mpi_spec,
        );
    else
        sim = @time MISMIP_PLUS();
    end 
	print(summary(sim))
end

