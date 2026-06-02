#add each of the individual ice thermodynamics models
include("./quadratic_temperature_approximation.jl")

struct NoThermoDynamics{T <: Real} <: AbstractThermoDynamics
    m :: T #uniform melt rate applied for grounded areas
end

NoThermoDynamics(; m = 0.0) = NoThermoDynamics(m)

"""
NoThermoDynamics(; <kwargs>)


Keyword arguments
=================
- basal_temperature  : basal ice temperature (K)
"""

"""
            update_ice_temperature_grounded_melt_rate!(thermo_dynamics::NoThermoDynamics, model::AbstractModel)

no thermodynamics model is used. ice temperature is constant in time
"""

function update_ice_temperature_grounded_melt_rate!(thermo_dynamics::NoThermoDynamics, model::AbstractModel)
    @unpack gh=model.fields
    gh.grounded_basal_melt[gh.mask] .= thermo_dynamics.m .* gh.grounded_fraction[gh.mask]
    return model
end