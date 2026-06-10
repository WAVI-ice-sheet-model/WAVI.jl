# WAVIBenchmarks — benchmark harness for WAVI.jl (BasicSpec / ThreadedSpec / MPISpec).
#
# To install dependencies the first time, or on Project.toml changes, run:
#   cd benchmarks
#   julia --project=. -e 'using Pkg; Pkg.instantiate()'
#
# Load in Julia:
#   using WAVIBenchmarks

module WAVIBenchmarks

using Comonicon

const BENCH_ROOT = dirname(@__DIR__)

# Include shared harness components
include(joinpath(BENCH_ROOT, "drivers.jl"))
include(joinpath(BENCH_ROOT, "utils.jl"))
include(joinpath(BENCH_ROOT, "plotting.jl"))
include(joinpath(@__DIR__, "harness.jl"))

export BenchmarkOptions, run_benchmark

"""
List registered benchmark driver adaptors discovered under `benchmarks/drivers/`.

This command checks the registry for available drivers that can be used 
with the benchmark harness and prints them to the console.
"""
@cast function drivers()
    names = available_drivers()
    println("Available drivers: ", isempty(names) ? "(none)" : join(names, ", "))
end

# Initialise Comonicon CLI
Comonicon.@main

# WIP: Run using `julia --project=benchmarks/ -e 'using WAVIBenchmarks; WAVIBenchmarks.drivers()`

end
