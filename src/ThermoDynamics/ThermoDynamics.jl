#add each of the individual ice thermodynamics models
include("./quadratic_temperature_approximation.jl")
include("./quadratic_temperature_approximation_iceberg_test.jl")

struct NoThermoDynamics <: AbstractThermoDynamics
end

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