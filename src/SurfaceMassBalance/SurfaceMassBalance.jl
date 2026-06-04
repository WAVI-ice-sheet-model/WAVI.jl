module SurfaceMassBalance

export update_accumulation_rate!

using WAVI: AbstractSurfaceMassBalance, AbstractModel, AbstractClimateForcing
using Parameters
using WAVI.Time: Clock
using WAVI.Grids

#return the sceme used to compute surface mass balance
get_surface_mass_balance(model::AbstractModel) = model.surface_mass_balance

#dispatch to specialised methods depending on surface mass balance law
"""

    update_accumulation_rate!(model::AbstractModel, clock::Clock)

Update the accumulation rate in the model using the time from clock.
"""
update_accumulation_rate!(model::AbstractModel, clock::Clock) = update_accumulation_rate!(get_surface_mass_balance(model)::AbstractSurfaceMassBalance, model, clock)

"""

    update_climate_forcing!(surface_mass_balance::AbstractSurfaceMassBalance) 

Generic wrapper function for updating the climate forcing. Overload this in your surface mass balance module of choice (see ISMIP7SMB.jl for an example)
"""
function update_climate_forcing!(surface_mass_balance::SMB, grid::Grid, clock::Clock) where {SMB <: AbstractSurfaceMassBalance}
    return nothing
end

function Grids.reconstruct_on_grid(smb::SMB, grid::Grid) where {SMB <: AbstractSurfaceMassBalance}
    return smb
end

function Grids.reconstruct_on_subdomain(smb::SMB, grid::Grid,subdomain::NTuple{4,<: Integer}) where {SMB <: AbstractSurfaceMassBalance}
    return smb
end

#include all specialised methods for computing accumulation rate
include("./AccumulationFromParams.jl")
include("./UniformAccumulationRate.jl")
include("./BinfileAccumulationRate.jl")
include("./ISMIP7SMB.jl")

#default behaviour if no clock passed
update_accumulation_rate!(model::AbstractModel) = update_accumulation_rate!(model::AbstractModel, Clock()) 

end