using WAVI

include("utils.jl")

function run_test()
    dx, dy = 8.e3, 8.e3;
    sx, sy = 640*1e3, 80*1e3;
    nx, ny = round(Int, sx/dx), round(Int, sy/dy);
    
    conditions = InitialConditions(
        initial_thickness = 20. .* ones(nx, ny),
        initial_u_veloc = [150 * sin(2x/nx) + 50 * tanh(1/ℯ^(1/y)) ^ 7 for x in 1:nx+1, y in 1:ny],
        initial_v_veloc = [30 - 5 * atan(x/nx) ^ 2 + 4 * cos((13/12(nx/x))/(11/0.5y)) ^ 7 for x in 1:nx, y in 1:ny+1],
    )

    @info "Preplotting conditions"
    plot_field(conditions.initial_u_veloc, "u.init.png")
    plot_field(conditions.initial_v_veloc, "v.init.png")

    u_iszero = falses(nx+1,ny);     # build x-direction velocity boundary condition matrix with no zero boundary conditions anywhere 
    u_iszero[1,:] .= true;          # set the x-direction velocity to zero at x = 0.
    
    v_iszero = falses(nx,ny+1);     # build x-direction velocity boundary condition matrix with no zero boundary conditions anywhere 
    v_iszero[:,1] .= true;          # set the y-direction velocity to zero at y = 0 (free slip)
    v_iszero[:,end] .= true;        # set the y-direction velocity to zero at y = 84km (free slip)
    v_iszero[1,:] .= true;          # set the y-direction velocity to zero at x = 0km (no slip in combination with u_iszero)

    grid = Grid(nx = nx, 
                ny = ny,   
                dx = dx, 
                dy = dy,
                u_iszero = u_iszero, 
                v_iszero = v_iszero)

    model = Model(grid, mismip_plus_bed, 100; initial_conditions=conditions, solver_params=SolverParams(maxiter_picard=1))

    plot_field(model.fields.gu.u, "u.before.png")
    plot_field(model.fields.gv.v, "v.before.png")
    update_state!(model)
    plot_field(model.fields.gu.u, "u.updated.png")
    plot_field(model.fields.gv.v, "v.updated.png")
    
    return model
end

run_test()
