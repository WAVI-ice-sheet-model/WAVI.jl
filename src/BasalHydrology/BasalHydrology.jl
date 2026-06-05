module BasalHydrology

export update_basal_water_thickness_effective_pressure!

using Parameters, Enzyme

using WAVI: AbstractBasalHydrology, AbstractModel
using WAVI.Grids


#add each of the individual basal hydrology models
include("./NoHydrology.jl")
include("./constant_basal_water_thickness.jl")
include("./leaky_bucket.jl")
include("./sheet_only_GlaDS.jl")

# Default behaviour. Overloaded for types that store spatial information.
function Grids.reconstruct_on_grid(basal_hydrology::BH, grid::Grid) where {BH <: AbstractBasalHydrology}
    return basal_hydrology
end

# Default behaviour. Overloaded for types that store spatial information.
function Grids.reconstruct_on_subdomain(basal_hydrology::BH, grid::Grid,subdomain::NTuple{4,<: Integer}) where {BH <: AbstractBasalHydrology}
    return basal_hydrology
end

end


