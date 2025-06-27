using WAVI
using MPI

ENV["GKSwstype"]="nul"

function mismip_plus_bed(x,y)
    xbar = 300000.0
    b0 = -150.0; b2 = -728.8; b4 = 343.91; b6 = -50.75;
    wc = 24000.0; fc = 4000.0; dc = 500.0;
    bx(x)= b0 + b2 * (x / xbar) ^ 2 + b4 * (x / xbar) ^ 4 + b6 * (x / xbar) ^ 6;
    by(y)= dc *( (1 + exp(-2(y - wc) / fc)) ^ (-1) + (1 + exp(2(y + wc) / fc)) ^ (-1) );
    b = max(bx(x) + by(y), -720.0);
    return b;
end

function plot_field(data, name="test.png")
    # Now we can plot
    Plots.heatmap(data);
    plot!(size = (600, 400), fmt=:png)
    savefig(name)
end

function run_test()
    dx, dy = 8.e3, 8.e3;
    sx, sy = 160.e3, 80.e3;
    nx, ny = round(Int, sx/dx), round(Int, sy/dy);

    u_iszero = falses(nx+1,ny);
    u_iszero[1,:] .= true;
    
    v_iszero = falses(nx,ny+1);
    v_iszero[:,1] .= true;
    v_iszero[:,end] .= true;
    v_iszero[1,:] .= true;

    grid = Grid(nx = nx, 
                ny = ny,   
                dx = dx, 
                dy = dy,
                u_iszero = u_iszero, 
                v_iszero = v_iszero)
    spec = MPISpec(1, 2, 3, grid)

    conditions = InitialConditions(
        initial_thickness = 20. .* ones(nx, ny),
        initial_u_veloc = [150 * sin(2x/nx) + 50 * tanh(1/ℯ^(1/y)) ^ 7 for x in 1:nx+1, y in 1:ny],
        initial_v_veloc = [30 - 5 * atan(x/nx) ^ 2 + 4 * cos((13/12(nx/x))/(11/0.5y)) ^ 7 for x in 1:nx, y in 1:ny+1],
    )

    # This creates the submodel from the global grid defined on a per-process basis
    model = Model(grid, mismip_plus_bed, 100, spec; 
        initial_conditions=conditions, solver_params=SolverParams(maxiter_picard=1))

    @info "Preplotting fields"
    u, v = model.global_fields.gu.u, model.global_fields.gv.v
    plot_field(model.fields.gu.u, "local.u.before.$(spec.rank).png")
    plot_field(model.fields.gv.v, "local.v.before.$(spec.rank).png")
    if spec.rank == 0
        plot_field(u, "global.u.before.$(spec.rank).png")
        plot_field(v, "global.v.before.$(spec.rank).png")
    end

    update_state!(model)

    @info "Postplotting fields"
    plot_field(model.fields.gu.u, "local.u.update.$(spec.rank).png")
    plot_field(model.fields.gv.v, "local.v.update.$(spec.rank).png")
    u, v = model.global_fields.gu.u, model.global_fields.gv.v
    if spec.rank == 0
        plot_field(u, "global.u.update.$(spec.rank).png")
        plot_field(v, "global.v.update.$(spec.rank).png")    
    end

    return model
end

run_test()
