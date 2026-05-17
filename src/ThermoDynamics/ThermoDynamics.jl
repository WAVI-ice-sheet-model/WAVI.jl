module ThermoDynamics

export update_ice_temperature_grounded_melt_rate!

using Parameters

using WAVI: AbstractThermoDynamics, AbstractModel

#add each of the individual ice thermodynamics models
include("./NoThermoDynamics.jl")
include("./quadratic_temperature_approximation.jl")

end