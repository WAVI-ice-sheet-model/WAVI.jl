module ThermoDynamics

import WAVI.Grids: reconstruct_on_grid, reconstruct_on_subdomain
import WAVI.ClimateForcing: update_climate_forcing!
export reconstruct_on_grid, reconstruct_on_subdomain,  update_climate_forcing!
export update_ice_temperature_and_basal_melt_rate!

using Parameters
using WAVI.Grids: Grid
using WAVI: AbstractThermoDynamics, AbstractModel
using WAVI.Grids
using WAVI.Time


#add each of the individual ice thermodynamics models
include("./NoThermoDynamics.jl")
include("./quadratic_temperature_approximation.jl")

end