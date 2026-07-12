# WAVIBenchmarks — benchmark harness for WAVI.jl (BasicSpec / ThreadedSpec / MPISpec).
#
# To install dependencies the first time, or on Project.toml changes, run:
#   cd benchmarks
#   julia --project=. -e 'using Pkg; Pkg.instantiate()'

module WAVIBenchmarks

using Comonicon

const BENCH_ROOT = dirname(@__DIR__)

# Include shared harness components
include(joinpath(BENCH_ROOT, "drivers.jl"))
include(joinpath(BENCH_ROOT, "resource_monitor.jl"))
include(joinpath(BENCH_ROOT, "utils.jl"))
include(joinpath(BENCH_ROOT, "plotting.jl"))
include(joinpath(@__DIR__, "harness.jl"))

export BenchmarkOptions, run_benchmark, set_benchmark_command!

"""
List registered benchmark driver adaptors discovered under `benchmarks/drivers/`.

This command checks the registry for available drivers that can be used 
with the benchmark harness and prints them to the console.
"""
@cast function drivers()
    names = available_drivers()
    println("Available drivers: ", isempty(names) ? "(none)" : join(names, ", "))
end

"""
Run a timed benchmark for the given execution mode and driver adaptor.

# Arguments

- `mode`: `basic`, `threaded`, or `mpi`
- `driver`: registered adaptor name (e.g. `mismip_plus`)

# Options

- `--niterations <n>`: [ThreadedSpec, MPISpec] Schwarz/PoU solver iterations (default: 2)

- `--ngridsx <n>`: [ThreadedSpec] x domain decomposition (default: 2)
- `--ngridsy <n>`: [ThreadedSpec] y domain decomposition (default: 2)
- `--overlap <n>`: [ThreadedSpec] Schwarz overlap cells (default: 2)

- `--px <n>`: [MPISpec] process grid x dimension (default: MPI world size)
- `--py <n>`: [MPISpec] process grid y dimension (default: 1)

- `--sample-interval <s>`: resource sample period in seconds (default: 0.25)
- `--no-plots`: skip NetCDF plots
- `--warmup`: run once untimed before the measured run
"""
@cast function run(
    mode::String,
    driver::String;
    ngridsx::Int = 2,
    ngridsy::Int = 2,
    overlap::Int = 2,
    niterations::Int = 2,
    px = nothing,
    py = nothing,
    sample_interval::Float64 = 0.25,
    no_plots::Bool = false,
    warmup::Bool = false,
)
    opts = BenchmarkOptions(
        mode;
        driver = driver,
        ngridsx = ngridsx,
        ngridsy = ngridsy,
        overlap = overlap,
        niterations = niterations,
        px = px isa Integer ? Int(px) : nothing,
        py = py isa Integer ? Int(py) : nothing,
        sample_interval = sample_interval,
        no_plots = no_plots,
        warmup = warmup,
    )
    run_benchmark(opts)
end

# Initialise Comonicon CLI
Comonicon.@main

end
