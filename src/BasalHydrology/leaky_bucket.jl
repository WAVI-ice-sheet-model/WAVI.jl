struct LeakyBucket{T <: Real} <: AbstractBasalHydrology
    drainage_rate :: T
    max_water_thickness :: T
    bed_roughness_scale :: T
    min_effective_pressure :: T
end

"""
LeakyBucket(; <kwargs>)


Keyword arguments
=================
- drainage_rate           : constant drainage rate (m/yr)
- max_water_thickness     : limits the basal water thickness (m)
- bed_roughness_scale     : effective bed roughness scale (m)
- min_effective_pressure  : minimum effective pressure (Pa)
"""

function LeakyBucket(;
                    drainage_rate = 0.01,
                    max_water_thickness = 10.,
                    bed_roughness_scale = 0.1,
                    min_effective_pressure = 1.0e4) 
                        
    return LeakyBucket(
                    drainage_rate,
                    max_water_thickness,
                    bed_roughness_scale,
                    min_effective_pressure)
end

"""
            update_basal_water_thickness_effective_pressure!(basal_hydrology::LeakyBucket, model::AbstractModel)

use a leaky bucket model to calculate the effective pressure
"""

function update_basal_water_thickness_effective_pressure!(basal_hydrology::LeakyBucket, model::AbstractModel)
    @unpack gh=model.fields
    @unpack params=model
    gh.basal_water_thickness .= max.(min.(gh.basal_water_thickness .+ (gh.grounded_basal_melt .- basal_hydrology.drainage_rate).*params.dt, basal_hydrology.max_water_thickness), 0.)
    
    # G. Flowers, 2000
    gh.effective_pressure .= max.(basal_hydrology.min_effective_pressure, params.density_ice .* params.g .* gh.h .* (1.0 .- min.(gh.basal_water_thickness ./ basal_hydrology.bed_roughness_scale, 1.0).^(3.5)))  
    gh.effective_pressure .= gh.effective_pressure .* gh.grounded_fraction
    
    return model
end