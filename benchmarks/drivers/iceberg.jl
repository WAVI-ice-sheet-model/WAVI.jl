# # Adaptor pointing to the Iceberg example driver.
# This script's filename (`iceberg.jl`) MUST match the module name below.

module iceberg

using MPI
using WAVI

# Include the main Iceberg example driver script
# Path is relative to this script
include(joinpath(@__DIR__, "..", "..", "example_drivers", "Iceberg", "Iceberg.jl"))

# # Required exports for the benchmark harness
const driver_run = Iceberg

# Function returning the global Grid
const driver_grid = Iceberg_Grid

# Vector of NetCDF variable names for optional post-run visualisation
const driver_plot_vars = ["h", "u", "v"]

end
