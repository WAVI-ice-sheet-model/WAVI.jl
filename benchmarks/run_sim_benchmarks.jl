using ArgParse
using BenchmarkTools
using LinearAlgebra
using Logging
using WAVI

function parse_cli()
    arg_settings = ArgParseSettings(exit_after_help = true)
    @add_arg_table! arg_settings begin
        "--depth"
            help = "Vertical depth on grid"
            default = 4
            arg_type = Int
        "--grid-cells"
            help = "Override the number of grid cells square"
            default = "100"
        "--cell-spacing"
            help = "Override the test grid cell spacing"
            default = "8000.0"            
        "--skip-dist-basic"
            help = "Skip the basic parallel spec distribution strategy"
            action = :store_true
        "--skip-dist-sharedmem"
            help = "Skip the shared memory spec distribution strategy"
            action = :store_true
        "--verbose"
            help = "Turn logging debug messages on"
            action = :store_true
        "threads"
            help = "Number of threads to run basic and shared memory with"
            required = true
    end
    return parse_args(ARGS, arg_settings)
end

args = parse_cli()

thread_nums = parse.(Int, split(args["threads"], ","))
num_cells = parse.(Int, split(args["grid-cells"], ","))
cell_spacing = parse.(Float64, split(args["cell-spacing"], ","))
depth = Int(args["depth"])
distributions = []

if ~args["skip-dist-basic"]
    push!(distributions, BasicParallelSpec())
end

if ~args["skip-dist-sharedmem"]
    push!(distributions, SharedMemorySpec(ngridsx=4, ngridsy=4, niterations=4))
end

const suite = BenchmarkGroup()
suite["runs"] = BenchmarkGroup()

for threads in thread_nums, cells in num_cells, spacing in cell_spacing, dist_spec in distributions
    suite["runs"][typeof(dist_spec)][threads][(cells, spacing, depth)] = @benchmarkable begin
        BLAS.set_num_threads($threads)
        min_elev = 20
        max_elev = 100
        
        @info "Creating test scenario with grid: $($cells) x $($cells) x $($depth) with $($spacing) spacing, $($threads) threads, distribution $(typeof($dist_spec))"
        
        grid = Grid(
            nx = $cells,
            ny = $cells,
            dx = $spacing,
            dy = $spacing,
            nσ = $depth,
            x0 = 0.0,
            y0 = -40000.0,
        )
        @info grid
        init = InitialConditions()
        params = Params()
        timestep = TimesteppingParams(
            dt = 1.0,
            end_time= 10., 
        )
        span = range(min_elev, max_elev, length = $cells)
            
        model = Model(
            grid = grid, 
            bed_elevation = [y - cos(ℯ/ (y/length(span))) + 2 * sin(x * 1e-3 * length(span)) for x in span, y in span],
            params = params,
            solver_params = SolverParams(),
            initial_conditions = init,
            melt_rate = UniformMeltRate(),
            parallel_spec = $dist_spec)

        sim = Simulation(
            model = model,
            output_params = OutputParams(),
            timestepping_params = timestep)

        run_simulation!(sim);        
    end
end

# If a cache of tuned parameters already exists, use it, otherwise, tune and cache
# the benchmark parameters. Reusing cached parameters is faster and more reliable
# than re-tuning `suite` every time the file is included.
paramspath = joinpath(dirname(@__FILE__), "params.json")

if isfile(paramspath)
    loadparams!(suite, BenchmarkTools.load(paramspath)[1], :evals)
else
    tune!(suite)
    BenchmarkTools.save(paramspath, params(suite))
end

results = run(suite, verbose = true)

@show results