```@meta
EditURL = "<unknown>/examples/overdeepened_bed.jl"
```

# MISMIP+ example

This example demonstrates the Marine Ice Sheet Model Intercomparison (MISMIP) + ice0 experiment (doi: 10.5194/tc-14-2283-2020).
This intercomparison exercise considers a plan-view ice sheet, in which the grounding line can stabilize on a section of bed which has a locally positive slope in the flow direction.
Such configurations are theoretically impossible in one horizontal dimension (see doi: 10.1029/2006JF000664), demonstrating the importance of variations in the second dimension for buttressing ice sheets.
This example demonstrates how to apply boundary conditions in `WAVI.jl`, and control the number of iterations in the velocity solve, and zip the output into a friendly format.

## Install dependencies

First let's make sure we have all required packages installed.

using Pkg
Pkg.add("https://github.com/RJArthern/WAVI.jl"), Pkg.add(Plots)

````@example overdeepened_bed
using WAVI, Plots, NCDatasets
````

## Basal Topography

The MISMIP+ domain is 640km in the x-direction and 80km in the y-direction, centred around $y = 0$.
The basal topography is given by $z_b = \max [B_x(x) + B_y(y), -720] where $B_x(x)$ is a sixth order, even polynomial and $B_y(y)$ introduces two bumps in the domain:

````@example overdeepened_bed
function mismip_plus_bed(x,y)
    xbar = 300000.0
    b0 = -150.0; b2 = -728.8; b4 = 343.91; b6 = -50.75;
    wc = 24000.0; fc = 4000.0; dc = 500.0;
    bx(x)=b0+b2*(x/xbar)^2+b4*(x/xbar)^4+b6*(x/xbar)^6;
    by(y)= dc*( (1+exp(-2(y-wc)/fc))^(-1) + (1+exp(2(y+wc)/fc))^(-1) );
    b = max(bx(x) + by(y), -720.0);
    return b;
end
````

Let's take a look at this bed. First we define the grid sizes and build some arrays, so we can plot. We'll use a high resolution to get a nice plot:

````@example overdeepened_bed
dx = 1.e3;
dy = 1.e3;
nx = round(Int, 640*1e3/dx);
ny = round(Int, 80*1e3/dx);
xx=[i*dx for i=1:nx, j=1:ny];
yy=[j*dy for i=1:nx, j=1:ny] .- 42000;
x = xx[:,1];
y = yy[1,:];
nothing #hide
````

Now we can plot

````@example overdeepened_bed
plt =  Plots.heatmap(x/1e3, y/1e3, mismip_plus_bed.(xx,yy)',
                    xlabel = "x (km)",
                    ylabel = "y (km)",
                    colorbar_title = "bed depth (m)")
plot!(size = (800,400))
````

## Boundary Conditions

The MISMIP+ experiment specifies no slip (zero velocity in both directions) boundary conditions at $x = 0$, and free-slip boundary conditions (zero velocity in the direction normal to the walls) at the lateral boundaries at $y = 0$km and $y = 84$km.
First, let's redefine the grid size to be lower resolution (to make the later simulations quicker)

````@example overdeepened_bed
dx = 4.e3;
dy = 4.e3;
nx = round(Int, 640*1e3/dx);
ny = round(Int, 80*1e3/dx);
nothing #hide
````

Velocity boundary conditions are controlled by specifying zeros in appropriate entries in arrays, which are then passed to the grid:

````@example overdeepened_bed
u_iszero = falses(nx+1,ny); #build x-direction velocity boundary condition matrix with no zero boundary conditions anywhere
u_iszero[1,:].=true;        #set the x-direction velocity to zero at x = 0.
v_iszero=falses(nx,ny+1);   #build x-direction velocity boundary condition matrix with no zero boundary conditions anywhere
v_iszero[:,1].=true;        #set the y-direction velocity to zero at y = 0 (free slip)
v_iszero[:,end].=true;       #set the y-direction velocity to zero at y = 84km (free slip)
v_iszero[1,:].=true;         #set the y-direction velocity to zero at x = 0km (no slip in combination with u_iszero)
nothing #hide
````

Now we build the grid as usual, passing the arrays we just constructed via optional arguments.

````@example overdeepened_bed
grid = Grid(nx = nx,
            ny = ny,
            dx = dx,
            dy = dy,
            u_iszero = u_iszero,
            v_iszero = v_iszero);
nothing #hide
````

## Solver Parameters
We're interested in the steady state reached at late times, rather than the solution along the way. We don't need to get the velocity right along the way, just have it correct eventually.
We therefore set the number of iterations in the velocity solve to be small: at each timestep, the solver just does a small number of iterations, and the velocity is only approximate. But, since we do a lot of iterations getting to steady state, the velocity gets to the right thing eventually.
This number, and other parameters relating to how the equations are solved, are set via a `SolverParams` object:

````@example overdeepened_bed
solver_params = SolverParams(maxiter_picard = 1);
nothing #hide
````

Explicitly, we set the number of iterations (formally, Picard iterations) to be as small as possible, i.e. one iteration.

## Make the model
Now we have our grid, bed, and solver parameters, we just need to set the appropriate initial conditions and physical parameters for MISMIP+, and then we can build our model.
In MISMIP+, the initial condition specifies 100m thick ice everywhere

````@example overdeepened_bed
initial_conditions = InitialConditions(initial_thickness = 100 .* ones(nx,ny));
nothing #hide
````

And the accumulation rate is set to 0.3m/a (all other defaults in WAVI are chosen according to the values in MISMIP+)

````@example overdeepened_bed
params = Params( accumulation_rate = 0.3);
nothing #hide
````

Now let's make our model! Note that we use the functional form of the bed (the array we plotted earlier has higher resolution than our model has)

````@example overdeepened_bed
model = Model(grid = grid,
            bed_elevation = mismip_plus_bed,
            initial_conditions = initial_conditions,
            params = params,
            solver_params = solver_params);
nothing #hide
````

## Assembling the Simulation
To get to steady state, we need to run our simulation for a good chunk of time, on the order of 10s of 1000s of years. We'll run for 10000 years, with a timestep of half a year:

````@example overdeepened_bed
timestepping_params = TimesteppingParams(dt = 0.5,
                                        end_time = 100.,);
nothing #hide
````

NB!! We have to specify the end time as `10000.` (a float number) rather than `10000` (an integer) because WAVI.jl expects the same type for the timestep `dt` and the end time `end_time`.

We'll output the solution along the way, and use this to convince ourselves later than we are in steady state. First let's make a directory to store the output

````@example overdeepened_bed
folder = joinpath(@__DIR__, "overdeepened_bed");
isdir(folder) && rm(folder, force = true, recursive = true);
mkdir(folder) ;
nothing #hide
````

Then define our output parameters. We'll output the thickness and grounded fraction every 200 years, and set the zip_format keyword argument to zip the output files.

````@example overdeepened_bed
output_params = OutputParams(outputs = (h = model.fields.gh.h,grfrac = model.fields.gh.grounded_fraction),
                            output_freq = 10.,
                            output_path = folder,
                            zip_format = "nc");
nothing #hide
````

Now we assemble our simulation, taking in the model, output parameters and timestepping_params:

````@example overdeepened_bed
simulation = Simulation(model = model, timestepping_params = timestepping_params, output_params = output_params);
nothing #hide
````

## Timestepping
Now all that's left to do is run our simulation! This is a long simulation and might take a while (~10 mins on my laptop)

````@example overdeepened_bed
run_simulation!(simulation);
nothing #hide
````

## Visualization
Let's have a look at the steady state thickness:

````@example overdeepened_bed
Plots.heatmap(simulation.model.grid.xxh[:,1]/1e3, simulation.model.grid.yyh[1,:]/1e3, simulation.model.fields.gh.h',
                xlabel = "x (km)",
                ylabel = "y (km)",
                colorbar_title = "ice thickness (m)")
````

And add the grounding line, which is where the grounded fraction transitions between 0 and 1 (grounded_fraction takes the value 1 at fully grounded grid points and 0 at fully floating grid points.) Our choice of 0.5 is somewhat arbitrary here -- any value between 0 and 1 will do!

````@example overdeepened_bed
Plots.contour!(simulation.model.grid.xxh[:,1]/1e3,
            simulation.model.grid.yyh[1,:]/1e3,
            simulation.model.fields.gh.grounded_fraction',
            fill = false,
            levels = [0.5,0.5],
            linecolor = :blue,
            linewidth = 2)
plot!(size = (1000,550))
````

You can see, by comparing with the plot of the bed earlier, that the grounding line sits on an overdeepened section of the bed!

Finally, let's check that it's in steady state, by looking at the evolution of the volume above floatation:

````@example overdeepened_bed
filename = joinpath(@__DIR__, "overdeepened_bed", "outfile.nc");

time, h, grfrac = NCDataset(filename, "r") do ds
    # Read variables directly from the NCDataset object
    h_data = ds["h"][:,:,:]
    grfrac_data = ds["grfrac"][:,:,:]
    time_data = ds["TIME"][:]
    return time_data, h_data, grfrac_data
end

vaf = sum(h .* grfrac .* dx .* dy, dims=(1, 2)); #volume of ice above floatation

Plots.plot(time, vaf[1,1,:]/1e9,
             marker = true,
             label = :none,
             xlabel = "time",
             ylabel = "volume above floatation (km^3)",
             framestyle = :box)
````

The volume above floatation reaches a plateau, suggesting that we have reached a steady state.
Finally, we clear up the files we just outputted

````@example overdeepened_bed
rm(folder, force = true, recursive = true);
nothing #hide
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

