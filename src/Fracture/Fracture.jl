module Fracture

export update_damage!, update_strain_history!

using Parameters

using WAVI: AbstractFracture, AbstractModel
using WAVI.Advection

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

include("./ConstantDamage.jl")
include("./DruckerPragerPhaseField.jl")
include("./utils.jl")

end