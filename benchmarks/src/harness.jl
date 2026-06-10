# # WAVI.jl benchmark harness

using MPI
using WAVI

const _VALID_MODES = (:basic, :threaded, :mpi)

# # Configuration

"""
    BenchmarkOptions

Configuration for a benchmark run.

# Fields

- `mode`: execution mode: `:basic`, `:threaded`, or `:mpi`
- `driver`: registered adaptor name (e.g. `"mismip_plus"`)
- `ngridsx`, `ngridsy`, `overlap`, `niterations`: ThreadedSpec / Schwarz parameters
- `px`, `py`: MPI process grid dimensions (`nothing` defaults to `Comm_size` × 1)
"""
Base.@kwdef struct BenchmarkOptions
    mode::Symbol
    driver::String
    ngridsx::Int = 2
    ngridsy::Int = 2
    overlap::Int = 2
    niterations::Int = 2
    px::Union{Int, Nothing} = nothing
    py::Union{Int, Nothing} = nothing
end

# Convenience constructor accepting a string mode (from the Comonicon CLI),
# coerced to the lower-case Symbol the dispatch expects.
BenchmarkOptions(mode::AbstractString; kwargs...) =
    BenchmarkOptions(; mode = Symbol(lowercase(mode)), kwargs...)

# Set from benchmarks/run.jl before the CLI runs; logged at the start of run_benchmark.
const BENCHMARK_COMMAND = Ref{Union{String, Nothing}}(nothing)

"""
    set_benchmark_command!(cmd::AbstractString)

Store the exact command-line invocation used to launch the benchmark.
This is logged to the console at the start of `run_benchmark` and is intended
to be saved alongside benchmark metadata for reproducibility.
"""
set_benchmark_command!(cmd::AbstractString) = (BENCHMARK_COMMAND[] = String(cmd))

# # Execution dispatch

"""
    run_benchmark(opts::BenchmarkOptions)

Load the requested driver adaptor, build appropriate WAVI spec, and run
`benchmark_main` for monitoring.
"""
function run_benchmark(opts::BenchmarkOptions)
    opts.mode in _VALID_MODES || error("Unknown mode '$(opts.mode)'. Use basic, threaded, or mpi.")

    # Initialise MPI up front so the entire setup is covered by finalisation.
    opts.mode == :mpi && !MPI.Initialized() && MPI.Init()

    rank = opts.mode == :mpi ? MPI.Comm_rank(MPI.COMM_WORLD) : 0
    rank == 0 && !isnothing(BENCHMARK_COMMAND[]) && @info "Command: $(BENCHMARK_COMMAND[])"

    driver = load_driver(opts.driver)

    try
        run_id = ""
        spec_kwargs = Dict{Symbol, Any}()

        if opts.mode == :basic
            # Serial: BasicSpec setup
            run_id = "basic"

        elseif opts.mode == :threaded
            # Shared memory: ThreadedSpec setup
            nx, ny, ov, ni = opts.ngridsx, opts.ngridsy, opts.overlap, opts.niterations
            spec_kwargs[:spec] = ThreadedSpec(; ngridsx = nx, ngridsy = ny, overlap = ov, niterations = ni)
            run_id = "threaded.$(nx)x$(ny)_o$(ov)_i$(ni)"

        elseif opts.mode == :mpi
            # Distributed: MPISpec setup
            comm = MPI.COMM_WORLD
            sz = MPI.Comm_size(comm)

            px = something(opts.px, sz)
            py = something(opts.py, 1)
            px * py == sz || error("MPI process grid px×py ($(px)×$(py)) must equal world size ($(sz)).")

            grid = Base.invokelatest(driver.grid)
            spec = MPISpec(px, py, 2, grid; pou = false, niterations = opts.niterations)

            spec_kwargs[:grid] = grid
            spec_kwargs[:spec] = spec
            run_id = "mpi.$(px)x$(py)_sz$(sz)"
        end

        benchmark_main(run_id, driver.run, spec_kwargs, driver.plot_vars, rank)
    finally
        opts.mode == :mpi && MPI.Initialized() && MPI.Finalize()
    end

    return nothing
end
