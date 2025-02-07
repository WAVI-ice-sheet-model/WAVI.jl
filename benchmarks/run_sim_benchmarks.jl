using BenchmarkTools
using LinearAlgebra
using WAVI

const suite = BenchmarkGroup()

suite["runs"] = BenchmarkGroup()

thread_nums = [1, 8, 16,]
grid_sizes = [(100, 8000.0)]#, (400, 2000.0), (1000, 800.0),]
depths = [4, 10]
distributions = [BasicParallelSpec(), SharedMemorySpec(ngridsx=4, ngridsy=4, niterations=4)]

for threads in thread_nums, grid_size in grid_sizes, depth in depths, dist_spec in distributions

    suite["runs"][typeof(dist_spec)][(threads, grid_size..., depth)] = @benchmarkable begin
        BLAS.set_num_threads($threads)
        min_elev = 20
        max_elev = 100
        n, d = $grid_size
        depth = $depth
        
        @info "Creating test scenario with grid: $n x $n x $depth with $d spacing, $($threads) threads, distribution $(typeof($dist_spec))"
        
        grid = Grid(
            nx = n,
            ny = n,
            dx = d,
            dy = d,
            nσ = depth,
            x0 = 0.0,
            y0 = -40000.0,
        )
        init = InitialConditions()
        params = Params()
        timestep = TimesteppingParams(
            dt = 1.0,
            end_time= 10., 
        )
        span = range(min_elev, max_elev, length = n)
            
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