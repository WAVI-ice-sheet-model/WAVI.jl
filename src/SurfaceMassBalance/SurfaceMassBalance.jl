module SurfaceMassBalance

export update_accumulation_rate!

using WAVI: AbstractSurfaceMassBalance, AbstractModel
using Parameters
using WAVI.Time: Clock

#return the sceme used to compute surface mass balance
get_surface_mass_balance(model::AbstractModel) = model.surface_mass_balance

#dispatch to specialised methods depending on surface mass balance law
"""

    update_accumulation_rate!(model::AbstractModel, clock::Clock)

Update the accumulation rate in the model using the time from clock.
"""
update_accumulation_rate!(model::AbstractModel, clock::Clock) = update_accumulation_rate!(get_surface_mass_balance(model)::AbstractSurfaceMassBalance, model, clock)

#include all specialised methods for computing accumulation rate
include("./AccumulationFromParams.jl")
include("./UniformAccumulationRate.jl")
include("./BinfileAccumulationRate.jl")

#default behaviour if no clock passed
update_accumulation_rate!(model::AbstractModel) = update_accumulation_rate!(model::AbstractModel, Clock()) 

end