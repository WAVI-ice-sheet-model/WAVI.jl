module ThermoDynamics

export update_ice_temperature_and_basal_melt_rate!

using Parameters

using WAVI: AbstractThermoDynamics, AbstractModel
using WAVI.Grids


#add each of the individual ice thermodynamics models
include("./NoThermoDynamics.jl")
include("./quadratic_temperature_approximation.jl")

# Default behaviour. Overloaded for types that store spatial information.
function Grids.reconstruct_on_grid(thermo_dynamics::TD, grid::Grid) where {TD <: AbstractThermoDynamics}
    return thermo_dynamics
end

# Default behaviour. Overloaded for types that store spatial information.
function Grids.reconstruct_on_subdomain(thermo_dynamics::TD, grid::Grid,subdomain::NTuple{4,<: Integer}) where {TD <: AbstractThermoDynamics}
    return thermo_dynamics
end

end