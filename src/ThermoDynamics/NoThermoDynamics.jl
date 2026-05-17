
export NoThermoDynamics

struct NoThermoDynamics <: AbstractThermoDynamics end

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
    return model
end