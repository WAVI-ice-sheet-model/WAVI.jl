"""
Iceberg(n_timesteps)

Driver routine for GlaDS verification test using WAVI.
Defines parameters, initialises state, and runs end_time_days/dt_days time steps.
"""
function GlaDS_test(;dt_days=0.1, end_time_days=5.0, melt_rate=7.93e-11, do_visu)

# GlaDS test geometry
nx, ny = 128, 4;
lx, ly = 100e3, 20e3;
dx, dy = lx / nx, ly / ny;
h_mask = trues(nx,ny);
hyd_potential_isfixed = falses(nx,ny);
hyd_potential_isfixed[1,:] .= true;

z_b = 0.0;  # bed elevation [m]
compute_H!(x) = 6.0 * (sqrt(x + 5000.0) - sqrt(5000.0)) + 1.0 - z_b; # Ice thickness function (1D, varies only with x)
xc = LinRange(0, lx, nx);
h = (zeros(nx, ny));
for j in 1:ny
    @. h[:, j] = compute_H!(xc);
end

bed_elevation = ones(nx,ny) .* z_b;

u_b = 1e-6;            # ice sliding speed       [m/s]
h_init = 0.05;         # water sheet thickness   [m]
basal_hydrology = SheetOnlyGlaDS(ncheck=10nx, maxiter=10nx^2, m=melt_rate);

dt = dt_days/365.25;
end_time = end_time_days/365.25;

# Set up boundary conditions
u_iszero = falses(nx+1,ny);  # build x-direction velocity boundary condition matrix with no zero boundary conditions anywhere
v_iszero = falses(nx,ny+1);  # build x-direction velocity boundary condition matrix with no zero boundary conditions anywhere
u_iszero[1,:]   .= true;     # set the x-direction velocity to zero at left boundary
u_iszero[end,:] .= true;     # set the x-direction velocity to zero at right boundary
v_iszero[:,1]   .= true;     # set the y-direction velocity to zero at bottom boundary
v_iszero[:,end] .= true;     # set the y-direction velocity to zero at top boundary

# Set up grid
grid = Grid(nx = nx,
            ny = ny,
            dx = dx,
            dy = dy,
            h_mask = h_mask,
            hyd_potential_isfixed = hyd_potential_isfixed,
            u_iszero = u_iszero,
            v_iszero = v_iszero);

initial_basal_water_thickness = h_init .* ones(nx, ny);
initial_conditions = InitialConditions(initial_thickness = h, initial_basal_water_thickness = initial_basal_water_thickness);

sliding_law = WeertmanSlidingLaw();
thermo_dynamics = NoThermoDynamics();
shelf_melt_rate = UniformMeltRate();
surface_mass_balance = AccumulationFromParams(),
params = Params(density_ice = 910.0, glen_a_ref=2.5e-25);
solver_params = SolverParams(maxiter_picard = 1);

model = Model(grid = grid,
            bed_elevation = bed_elevation,
            initial_conditions = initial_conditions,
            solver_params = SolverParams(maxiter_picard=30),
            params = params,
            basal_hydrology=basal_hydrology,
            sliding_law = sliding_law,
            shelf_melt_rate = shelf_melt_rate,
            surface_mass_balance = surface_mass_balance,
            thermo_dynamics = thermo_dynamics);

timestepping_params = TimesteppingParams(dt = dt, end_time = end_time);
simulation = Simulation(model = model, timestepping_params = timestepping_params);

t = 0;
while round(t*365.25, digits=2) < end_time_days
    model.fields.gh.bed_speed.=u_b.*model.params.sec_per_year; # [m/yr]
    model.fields.gh.h.=h;
    timestep!(simulation);
    t = simulation.clock.time;
    # println("t=",round(t*365.25, digits=2),", end_time=",round(end_time*365.25, digits=2)," [days]")
end

# plot
if do_visu
    @views av(x) = @. 0.5 * (x[1:end-1] + x[2:end])
    xv, yv = LinRange(dx / 2, lx - dx / 2, nx + 1), LinRange(dy / 2, ly - dy / 2, ny + 1)
    xc, yc = av(xv), av(yv)

    fig = Figure(; size=(1000, 500))

    # Choose slice indices for plotting
    y_slice_idx = Int(ny ÷ 2);  # middle slice in y-direction
    xc_slice = LinRange(dx / 2, lx - dx / 2, nx);

    # Create axis for ice thickness H (slice through middle in y-direction)
    ax_H = Axis(fig[1, 1];
                xlabel="Distance x [m]",
                ylabel="Ice Thickness H [m]",
                title="Ice Thickness (y-slice at y=$(yc[y_slice_idx]))")
    lines!(ax_H, xc_slice, simulation.model.fields.gh.h[:, y_slice_idx]; linewidth=2, color=:blue, label="Ice Thickness")

    # Create axis for effective pressure (slice through middle in y-direction)
    ax_N = Axis(fig[2, 1];
                xlabel="Distance x [m]",
                ylabel="Effective Pressure N [Pa]",
                title="Effective Pressure (t = $(round(t*365.25, digits=2)) days)")
    lines!(ax_N, xc_slice, simulation.model.fields.gh.effective_pressure[:, y_slice_idx]; linewidth=2, color=:green, label="N")

    # Create axis for hydraulic potential (slice through middle in y-direction)
    ax_ϕ = Axis(fig[1, 2];
                xlabel="Distance x [m]",
                ylabel="Hydraulic Potential ϕ [Pa]",
                title="Hydraulic Potential (t = $(round(t*365.25, digits=2)) days)")
    lines!(ax_ϕ, xc_slice, simulation.model.fields.gh.hydraulic_potential_b[:, y_slice_idx]; linewidth=2, color=:orange, label="ϕ")

    # Create axis for water sheet thickness (slice through middle in y-direction)
    ax_h = Axis(fig[2, 2];
                xlabel="Distance x [m]",
                ylabel="Water Sheet Thickness h [m]",
                title="Water Sheet Thickness (t = $(round(t*365.25, digits=2)) days)")
    lines!(ax_h, xc_slice, simulation.model.fields.gh.basal_water_thickness[:, y_slice_idx]; linewidth=2, color=:cornflowerblue, label="h")

    # Get data arrays at melting parameter m
    include("verification_tests/check-vs-gladsog.jl")
    ref_phi, ref_h = dataOG[melt_rate]

    # Find the closest reference time point
    if isapprox(t*365.25, 0; atol=dt_days)
        t_idx = 1;
    elseif isapprox(t*365.25, 5; atol=dt_days)
        t_idx = 2;
    elseif isapprox(t*365.25, 1000; atol=dt_days)
        t_idx = 3;
    else
        # Find the closest reference time point
        ref_times = [0, 5, 1000];
        _, t_idx = findmin(abs.(ref_times .- t*365.25));
        @warn "No exact match for time t=$(t*365.25) days. Using closest reference time: $(ref_times[t_idx]) days."
    end

    # Styling
    ref_colors  = [:darkblue, :darkred, :darkgreen];
    markers     = [:circle, :utriangle, :diamond];
    time_labels = ["Initial (t=0)", "Transient (t=5 days)", "Steady (t=1000 days)"];

    out_x       = 1e3 * [5, 15, 25, 35, 45, 55, 65, 75, 85, 95];
    scatter!(ax_ϕ, out_x, ref_phi[:, t_idx];
                color=ref_colors[t_idx], marker=markers[t_idx], markersize=12,
                strokewidth=1, strokecolor=:black,
                label="Ref: " * time_labels[t_idx])

    scatter!(ax_h, out_x, ref_h[:, t_idx];
                color=ref_colors[t_idx], marker=markers[t_idx], markersize=12,
                strokewidth=1, strokecolor=:black,
                label="Ref: " * time_labels[t_idx])

    # Add legends
    axislegend(ax_H; position=:rb)
    axislegend(ax_N; position=:rb)
    axislegend(ax_ϕ; position=:rb)
    axislegend(ax_h; position=:lb)

    # Update layout and display
    fig[0, :] = Label(fig, "SHMIP Suite A (2D C-grid, y-slice)\nm = $melt_rate"; fontsize=20)
    plotPath = "verification_tests/WAVI_GlaDS_test.png"
    println("save figure at ",plotPath)
    save(plotPath, fig)
end

return simulation

end
