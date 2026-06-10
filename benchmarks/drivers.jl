# # Registry and loader for benchmark driver adaptors
#
# This module registers and loads benchmark driver adaptors from
# `benchmarks/drivers/<name>.jl`. Each adaptor is a module
# (e.g. `mismip_plus`) inside `BenchmarkDrivers`. Adaptors wrap
# example drivers behind a common interface so the benchmarking
# harness is driver agnostic.
#
# Usage (from main WAVI.jl directory):
#   julia --project=. benchmarks/drivers.jl
#   julia --project=. -e 'include("benchmarks/drivers.jl"); load_driver("mismip_plus")'

module BenchmarkDrivers
end

const DRIVERS_DIR = joinpath(@__DIR__, "drivers")

# Cache for loaded drivers to avoid redefining modules
const _LOADED_DRIVERS = Dict{String, NamedTuple}()

"""
    available_drivers() -> Vector{String}

Return sorted driver names (from `benchmarks/drivers/<name>.jl` as `<name>`).
The module <name> inside each file must match the filename as well.
"""
function available_drivers()
    names = String[]
    for file in readdir(DRIVERS_DIR)
        endswith(file, ".jl") || continue
        push!(names, replace(file, ".jl" => ""))
    end
    return sort(names)
end

"""
    load_driver(name::String) -> NamedTuple

Load `benchmarks/drivers/<name>.jl` and return a named tuple containing
the `driver_run`, `driver_grid`, and `driver_plot_vars` functions.
"""
function load_driver(name::String)
    # Return from cache if already loaded
    if haskey(_LOADED_DRIVERS, name)
        return _LOADED_DRIVERS[name]
    end

    # Ensure requested driver file exists
    path = joinpath(DRIVERS_DIR, "$(name).jl")
    isfile(path) || error(
        "Unknown benchmark driver '$(name)'. " *
        "Available: $(join(available_drivers(), ", "))",
    )

    # Include driver adaptor into BenchmarkDrivers module
    mod_sym = Symbol(name)
    if !isdefined(BenchmarkDrivers, mod_sym)
        Base.include(BenchmarkDrivers, path)
    end

    # Use invokelatest because driver is defined dynamically at runtime
    result = invokelatest() do
        driver_mod = getfield(BenchmarkDrivers, mod_sym)

        # Verify adaptor contract is met
        for sym in (:driver_run, :driver_grid, :driver_plot_vars)
            isdefined(driver_mod, sym) || error(
                "Driver module \$(mod_sym) hasn't defined \$$(sym)",
            )
        end

        # Extract fields and return as NamedTuple
        return (
            run = getfield(driver_mod, :driver_run),
            grid = getfield(driver_mod, :driver_grid),
            plot_vars = getfield(driver_mod, :driver_plot_vars),
        )
    end

    # Store in cache
    _LOADED_DRIVERS[name] = result

    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    driver = load_driver("mismip_plus")
    println("Available drivers: ", join(available_drivers(), ", "))
    println("Loaded driver: mismip_plus")
    println("  run:       ", driver.run)
    println("  grid:      ", driver.grid)
    println("  plot_vars: ", driver.plot_vars)
end
