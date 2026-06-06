module ClimateForcing

export update_climate_forcing!

using WAVI: AbstractClimateForcing
using WAVI.Grids
using WAVI.Time

include("./ISMIP7.jl")

"""

    update_climate_forcing!(shelf_melt_rate::AbstractMeltRate, grid::Grid, clock::Clock) 

Generic wrapper function for updating the climate forcing. Overload this if needed.)
"""
function update_climate_forcing!(x, grid::Grid, clock::Clock)
    return nothing
end


end