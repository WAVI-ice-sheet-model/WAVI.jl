#add each of the individual basal hydrology models
include("./leaky_bucket.jl")

struct NoHydrology <: AbstractBasalHydrology
end

"""
            update_basal_water_thickness_effective_pressure!(basal_hydrology::NoHydrology, model::AbstractModel)

no basal hydrology model is used. effective pressure is constant in time
"""

function update_basal_water_thickness_effective_pressure!(basal_hydrology::NoHydrology, model::AbstractModel)
    @unpack gh=model.fields
    gh.effective_pressure .= model.params.effective_pressure .* gh.grounded_fraction
    return model
end