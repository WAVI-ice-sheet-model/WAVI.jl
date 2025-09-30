using WAVI 
function MISMIP_PLUS()
    #Grid and boundary conditions
    nx = 80
    ny = 10
    nσ = 4
    x0 = 0.0
    y0 = -40000.0
    dx = 8000.0
    dy = 8000.0
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
    solver_params = SolverParams(maxiter_picard = maxiter_picard)

    #Physical parameters
    default_thickness = 100.0 #set the initial condition this way
    accumulation_rate = 0.3
    default_temperature=265.700709
    params = Params(default_thickness = default_thickness, 
                    accumulation_rate = accumulation_rate,
                    default_temperature = default_temperature)

    #make the model
    model = Model(grid = grid,
                     bed_elevation = bed, 
                     params = params, 
                     solver_params = solver_params)

    #timestepping parameters
    niter0 = 0
    dt = 0.5
    end_time = 400.
    chkpt_freq = 1000.
    pchkpt_freq = 20.
    timestepping_params = TimesteppingParams(niter0 = niter0, 
                                            dt = dt, 
                                            end_time = end_time, 
                                            chkpt_freq = chkpt_freq, 
                                            pchkpt_freq = pchkpt_freq)

    #output parameters
    folder = "outputs"
    outputs = (h = model.fields.gh.h,
               u = model.fields.gh.u,
               v = model.fields.gh.v,
               b = model.fields.gh.b,
               grfrac = model.fields.gh.grounded_fraction) #output velocities and thickness
    output_freq = 20.
    output_params = OutputParams(outputs = outputs, 
                            output_path = folder,
                            output_freq = output_freq,
                            output_format = "jld2",
                            zip_format = "nc")
    
    simulation = Simulation(model = model, 
                        timestepping_params = timestepping_params,
                        output_params = output_params)
            
    #perform the simulation
    run_simulation!(simulation)
    return simulation
end

@time simulation = MISMIP_PLUS();