# # Adaptor pointing to the MISMIP+ example driver.
# This script's filename (`mismip_plus.jl`) MUST match the module name below.
#
# Each adaptor defines the following fields:
#   driver_run       - (; folder, grid, spec) -> simulation
#   driver_grid      - () -> global Grid (required for MPISpec)
#   driver_plot_vars - Vector of NetCDF variable names for optional post-run plots

module mismip_plus

using MPI
using WAVI

# Include the main MISMIP+ example driver script
# Path is relative to this script
include(joinpath(@__DIR__, "..", "..", "example_drivers", "MISMIP_PLUS", "MISMIP_PLUS.jl"))

# # Required exports for the benchmark harness
# Function returning a Simulation object
const driver_run = MISMIP_PLUS

# Function returning the global Grid (required for MPISpec initialisation)
const driver_grid = MISMIP_PLUS_GRID

# Vector of NetCDF variable names for optional post-run visualisation
const driver_plot_vars = ["h", "u", "v"]

end
