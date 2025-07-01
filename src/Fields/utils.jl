"""
    get_u_mask(h_mask)

Find mask of valid grid points on u-grid corresponding to a mask defined on h-grid.

"""
function get_u_mask(h_mask)
    #include all u faces next to a selected center
    (nx,ny)=size(h_mask)
    u_mask=falses(nx+1,ny)
    u_mask[1:end-1,1:end]=u_mask[1:end-1,1:end].|h_mask
    u_mask[2:end,1:end]=u_mask[2:end,1:end].|h_mask
    return u_mask
end
"""
    get_v_mask(h_mask)

Find mask of valid grid points on v-grid corresponding to a mask defined on h-grid.

"""
function get_v_mask(h_mask)
    #include all v faces next to a selected center
    (nx,ny)=size(h_mask)
    v_mask=falses(nx,ny+1)
    v_mask[1:end,1:end-1]=v_mask[1:end,1:end-1].|h_mask
    v_mask[1:end,2:end]=v_mask[1:end,2:end].|h_mask
    return v_mask
end
"""
    get_c_mask(h_mask)

Find mask of valid grid points on c-grid corresponding to a mask defined on h-grid.

"""
function get_c_mask(h_mask)
    #select cell corners with four neighbouring cell centers in h_mask
    c_mask=h_mask[1:end-1,1:end-1] .& h_mask[1:end-1,2:end] .& h_mask[2:end,1:end-1] .& h_mask[2:end,2:end]
    return c_mask
end

#1D Matrix operator utility functions.
spI(n) = spdiagm(n,n, 0 => ones(n))
∂1d(n,dx) = spdiagm(n,n+1,0 => -ones(n), 1 => ones(n))/dx
c(n) = spdiagm(n,n+1,0 => ones(n), 1 => ones(n))/2
χ(n) = spdiagm(n,n+2, 1 => ones(n))


"""
check_initial_conditions(initial_conditions, params)

Check whether initial conditions have been specified. Default them to standard values if not

TODO: should be encapsulated via outer constructors in some way
"""
function check_initial_conditions(initial_conditions, params, grid)
    @info "Setting up initial conditions based on params"
    if all(isnan.(initial_conditions.initial_thickness))
        default_thickness = params.default_thickness
        #@info "Did not find a specified initial thickness, reverting to default value specified in params ($default_thickness m everywhere)...\n...If you have set niter0 > 0 without invoking the update flag, you can ignore this message"
        initial_conditions = @set initial_conditions.initial_thickness =  default_thickness*ones(grid.nx, grid.ny)
    end

    if all(isnan.(initial_conditions.initial_grounded_fraction))
        #@info "Did not find a specified grounded fraction, reverting to default value of one everywhere)...\n...If you have set niter0 > 0 without invoking the update flag, you can ignore this message"
        initial_conditions = @set initial_conditions.initial_grounded_fraction =  ones(grid.nx, grid.ny)
    end

    if all(isnan.(initial_conditions.initial_u_veloc))
        #@info "Did not find a specified initial u velocity, reverting to default value of zero everywhere)...\n...If you have set niter0 > 0 without invoking the update flag, you can ignore this message"
        initial_conditions = @set initial_conditions.initial_u_veloc =  zeros(grid.nx+1, grid.ny)
    end

    if all(isnan.(initial_conditions.initial_v_veloc))
        #@info "Did not find a specified initial u velocity, reverting to default value of zero everywhere)...\n...If you have set niter0 > 0 without invoking the update flag, you can ignore this message"
        initial_conditions = @set initial_conditions.initial_v_veloc =  zeros(grid.nx, grid.ny+1)
    end

    if all(isnan.(initial_conditions.initial_viscosity))
        default_viscosity = params.default_viscosity
        #@info "Did not find a specified initial viscosity, reverting to default value specified in params ($default_viscosity Pa s everywhere)...\n...If you have set niter0 > 0 without invoking the update flag, you can ignore this message"
        initial_conditions = @set initial_conditions.initial_viscosity =  default_viscosity*ones(grid.nx, grid.ny, grid.nσ)
    end

    if all(isnan.(initial_conditions.initial_temperature))
        default_temperature = params.default_temperature
        #@info "Did not find a specified initial temperature, reverting to default value specified in params ($default_temperature K everywhere)...\n...If you have set niter0 > 0 without invoking the update flag, you can ignore this message"
        initial_conditions = @set initial_conditions.initial_temperature =  default_temperature*ones(grid.nx, grid.ny, grid.nσ)
    end

    if all(isnan.(initial_conditions.initial_damage))
        default_damage = params.default_damage
        #@info "Did not find a specified initial damage field, reverting to default value specified in params ($default_damage everywhere)...\n...If you have set niter0 > 0 without invoking the update flag, you can ignore this message"
        initial_conditions = @set initial_conditions.initial_damage =  default_damage*ones(grid.nx, grid.ny, grid.nσ)
    end



    #check sizes are compatible
    (size(initial_conditions.initial_thickness) == (grid.nx, grid.ny)) || throw(DimensionMismatch("Initial thickness field is not compatible with grid size. Input thickess field is has size $(size(initial_conditions.initial_thickness)), which must match horizontal grid size ($(grid.nx) x $(grid.ny))"))
    (size(initial_conditions.initial_grounded_fraction) == (grid.nx, grid.ny)) || throw(DimensionMismatch("Initial grounded fraction field is not compatible with grid size. Input grounded fraction field is has size $(size(initial_conditions.initial_grounded_fraction)), which must match horizontal grid size ($(grid.nx) x $(grid.ny))"))
    (size(initial_conditions.initial_u_veloc) == (grid.nx+1, grid.ny)) || throw(DimensionMismatch("Initial u_velocity field is not compatible with grid size. Input u_velocity field has size $(size(initial_conditions.initial_u_veloc)), which must match horizontal grid size ($(grid.nx+1) x $(grid.ny))"))
    (size(initial_conditions.initial_v_veloc) == (grid.nx, grid.ny+1)) || throw(DimensionMismatch("Initial v_velocity field is not compatible with grid size. Input v_velocity field has size $(size(initial_conditions.initial_v_veloc)), which must match horizontal grid size ($(grid.nx) x $(grid.ny+1))"))
    (size(initial_conditions.initial_temperature) == (grid.nx, grid.ny, grid.nσ)) || throw(DimensionMismatch("Initial temperature field is not compatible with grid size. Input temperature field is has size $(size(initial_conditions.initial_temperature)), which must match 3D grid size ($(grid.nx), $(grid.ny), $(grid.nσ))"))
    (size(initial_conditions.initial_viscosity) == (grid.nx, grid.ny, grid.nσ)) || throw(DimensionMismatch("Initial viscosity field is not compatible with grid size. Input temperature field is has size $(size(initial_conditions.initial_temperature)), which must match 3D grid size ($(grid.nx), $(grid.ny), $(grid.nσ))"))
    (size(initial_conditions.initial_damage) == (grid.nx, grid.ny, grid.nσ)) || throw(DimensionMismatch("Initial damage field is not compatible with grid size.Input temperature field is has size $(size(initial_conditions.initial_temperature)), which must match 3D grid size ($(grid.nx), $(grid.ny), $(grid.nσ))"))


    return initial_conditions
end
