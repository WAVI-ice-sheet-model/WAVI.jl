module BasalHydrology

export update_basal_water_thickness_effective_pressure!

using Parameters, Enzyme

using WAVI: AbstractBasalHydrology, AbstractModel


#add each of the individual basal hydrology models
include("./NoHydrology.jl")
include("./constant_basal_water_thickness.jl")
include("./leaky_bucket.jl")
include("./sheet_only_GlaDS.jl")

end


