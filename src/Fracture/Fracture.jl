module Fracture

export update_damage!, update_strain_history!

using Parameters
using NCDatasets

using WAVI: AbstractFracture, AbstractModel, AbstractGrid
using WAVI.Advection
using WAVI.Time: Clock
using WAVI.Grids: Grid



function reconstruct_on_grid(fracture::F, grid::Grid) where {F <: AbstractFracture}
    return fracture
end

function reconstruct_on_subdomain(fracture::F, grid::Grid,subdomain::NTuple{4,<: Integer}) where {F <: AbstractFracture}
    return fracture
end

get_fracture(model::AbstractModel{T,N}) where {T,N} = model.fracture

"""
    update_damage!(model::AbstractModel)

Update the damage field and the strain history from 'ice shelf' parts of strain rate tensor, neglecting all vertical shear.
"""
update_damage!(model::AbstractModel{T,N};kwargs ...) where {T,N} = update_damage!(get_fracture(model),model; kwargs ...)


"""
    update_strain_history!(model::AbstractModel)

Update the strain history stored in the model structure.
"""
update_strain_history!(model::AbstractModel{T,N};kwargs ...) where {T,N} = update_strain_history!(get_fracture(model),model; kwargs ...)


"""

    update_climate_forcing!(fracture::AbstractFracture,clock::Clock) 

Generic wrapper function for updating the climate forcing. Overload this in your fracture module of choice (see ISMIP7hydrofracture.jl for an example)
"""
function update_climate_forcing!(fracture::AbstractFracture, grid::AbstractGrid, clock::Clock)
    return nothing
end

include("./ConstantDamage.jl")
include("./DruckerPragerPhaseField.jl")
include("./ISMIP7hydrofracture.jl")
include("./utils.jl")

end