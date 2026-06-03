export ConstantBasalWaterThickness

struct ConstantBasalWaterThickness{T <: Real} <: AbstractBasalHydrology
    min_effective_pressure :: T
end

"""
ConstantBasalWaterThickness(; <kwargs>)


Keyword arguments
=================
- min_effective_pressure  : minimum effective pressure (Pa)
"""

function ConstantBasalWaterThickness(;
                    min_effective_pressure = 1.0e4) 
                        
    return ConstantBasalWaterThickness(
                    min_effective_pressure)
end

"""
            update_basal_water_thickness_effective_pressure!(basal_hydrology::ConstantBasalWaterThickness, model::AbstractModel; update_basal_water_thickness::Bool = true)

use a constant basal water thickness to calculate the effective pressure (assumed hydrostatic, so proportional to basal_water_thickness)
"""

function update_basal_water_thickness_effective_pressure!(basal_hydrology::ConstantBasalWaterThickness, model::AbstractModel; update_basal_water_thickness::Bool = true)
    @unpack gh=model.fields
    @unpack params=model
    gh.effective_pressure .= max.(basal_hydrology.min_effective_pressure, (params.density_ice .* params.g .* gh.h) .- (params.density_freshwater .* params.g .* gh.basal_water_thickness))
    gh.effective_pressure .= gh.effective_pressure .* gh.grounded_fraction
    return model
end