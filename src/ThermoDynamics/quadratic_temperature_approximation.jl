struct QuadraticTemperatureApproximation{T <: Real} <: AbstractThermoDynamics
    specific_heat_capacity :: T
    thermal_conductivity :: T
    latent_heat_fusion :: T
    geothermal_heat_flux :: T
    internal_heat_production :: T
    melting_point_coeff :: T
    kelvin_conversion :: T
    reg_thickness :: T
end

"""
QuadraticTemperatureApproximation(; <kwargs>)


Keyword arguments
=================
- specific_heat_capacity    : specific heat capacity of ice (J/(kg*K))
- thermal_conductivity      : thermal conductivity of ice (W/(m*K))
- latent_heat_fusion        : latent heat of fusion of ice (J/kg)
- geothermal_heat_flux      : geothermal heat flux (W/m²)
- internal_heat_production  : internal heat production (W/m³)
- melting_point_coeff       : dependence of melting point on pressure (K/Pa)
- kelvin_conversion         : convert kelvin to degreee celsius (K)
- reg_thickness             : regularization thickness, used to prevent ice thickness going to zero (m)
"""

function QuadraticTemperatureApproximation(;
                    specific_heat_capacity = 2.1e3, # ~2.1 J/(kg K) at 0°C (according to equation in Tarasov et al., 2025)
                    thermal_conductivity = 2.1, # ~2.1 W/(m K) or 66 MJ/(m yr K) at 0°C; ~2.6 W/(m K) or 83 MJ/(m yr K) at -60°C (Hooke, Principles of Glacier Mechanics, 3rd Edition)
                    latent_heat_fusion = 3.34e5, # J/kg (Hooke, Principles of Glacier Mechanics, 3rd Edition)
                    geothermal_heat_flux = 0.05, # world average (Hooke, Principles of Glacier Mechanics, 3rd Edition)
                    internal_heat_production = 0.,
                    melting_point_coeff = 9.76e-8, # (Rutt et al., 2009)
                    kelvin_conversion = 273.15,
                    reg_thickness=1.0e-5)
                        
    return QuadraticTemperatureApproximation(
                    specific_heat_capacity,
                    thermal_conductivity,
                    latent_heat_fusion,
                    geothermal_heat_flux,
                    internal_heat_production,
                    melting_point_coeff,
                    kelvin_conversion,
                    reg_thickness)
end

"""
            update_ice_temperature_grounded_melt_rate!(thermo_dynamics::QuadraticTemperatureApproximation, model::AbstractModel)

use a quadratic temperature approximation to calculate the ice temperature and grounded basal melt rate
also updates the stiffness parameter B in Glen flow law.
"""

function update_ice_temperature_grounded_melt_rate!(thermo_dynamics::QuadraticTemperatureApproximation{T}, model::AbstractModel) where {T}
    @unpack gh,gu,gv,g3d=model.fields
    @unpack params=model
    thermal_diffusivity = thermo_dynamics.thermal_conductivity / (params.density_ice * thermo_dynamics.specific_heat_capacity)

    # calculate depth-averaged temperature on u and v grid
    θ_ave_gu = zeros(T,gu.nxu,gu.nyu)
    θ_ave_gv = zeros(T,gv.nxv,gv.nyv)
    onesvec = ones(T,gh.nxh*gh.nyh)
    θ_ave_gu[gu.mask] .= (gu.samp*(gu.centᵀ*(gh.crop*gh.θ_ave[:])))./(gu.samp*(gu.centᵀ*(gh.crop*onesvec)))
    θ_ave_gv[gv.mask] .= (gv.samp*(gv.centᵀ*(gh.crop*gh.θ_ave[:])))./(gv.samp*(gv.centᵀ*(gh.crop*onesvec)))

    # calculate horizontal heat flux
    horizontal_heat_flux = zeros(T,gh.nxh,gh.nyh)
    horizontal_heat_flux[gh.mask] = gh.samp*((gu.∂x*(gu.crop*(gu.h[:].*gu.u[:].*θ_ave_gu[:]))) .+ 
                                    (gv.∂y*(gv.crop*(gv.h[:].*gv.v[:].*θ_ave_gv[:]))))./params.sec_per_year # K*m/s

    for j=1:g3d.nys
        for i=1:g3d.nxs
            if gh.mask[i,j]
                # calculate temperature balance terms
                θ_pressure_melting_base = thermo_dynamics.kelvin_conversion - params.density_ice * params.g * thermo_dynamics.melting_point_coeff * gh.h[i,j]
                θ_1_trial = -gh.h[i,j] / thermo_dynamics.thermal_conductivity * (thermo_dynamics.geothermal_heat_flux + gh.τbed[i,j] * (gh.bed_speed[i,j]/params.sec_per_year))
                θ_2 = g3d.θ[i,j,end] - g3d.θ[i,j,1] - θ_1_trial

                vertical_heat_flux = thermal_diffusivity * (2. * θ_2 / gh.h[i,j])
                internal_heat_source = thermo_dynamics.internal_heat_production * gh.h[i,j] / (params.density_ice * thermo_dynamics.specific_heat_capacity)

                # calculate depth-averaged temperature change dθ_ave_dt
                dθ_ave_dt = (gh.accumulation[i,j]/params.sec_per_year * g3d.θ[i,j,end] - 
                            gh.basal_melt[i,j]/params.sec_per_year * g3d.θ[i,j,1] + 
                            vertical_heat_flux + internal_heat_source -
                            gh.θ_ave[i,j] * gh.dhdt[i,j]/params.sec_per_year -
                            horizontal_heat_flux[i,j]) / gh.h[i,j]

                # update depth-averaged temperature θ_ave
                gh.θ_ave[i,j] = gh.θ_ave[i,j] + params.dt*params.sec_per_year * dθ_ave_dt
                
                # limit basal temperature to the pressure melting point
                θ_base_trial = (3/2 * gh.θ_ave[i,j]) - (1/4 * θ_1_trial) - (1/2 * g3d.θ[i,j,end])
                g3d.θ[i,j,1] = min(θ_pressure_melting_base,θ_base_trial)

                # adjust θ_1 and θ_2 to be in line with the basal temperature limited to the pressure melting point
                θ_1 = (6 * gh.θ_ave[i,j]) - (2 * g3d.θ[i,j,end]) - (4 * g3d.θ[i,j,1])
                θ_2 = g3d.θ[i,j,end] - g3d.θ[i,j,1] - θ_1

                # calculate grounded basal melt rate
                q_melt = -thermo_dynamics.thermal_conductivity / gh.h[i,j] * (θ_1_trial - θ_1) # W/m²
                gh.grounded_basal_melt[i,j] = q_melt / (params.density_ice * thermo_dynamics.latent_heat_fusion) * params.sec_per_year # m/yr
                gh.grounded_basal_melt[i,j] = gh.grounded_basal_melt[i,j] * gh.grounded_fraction[i,j]
                
                # calculate ice temperature profile
                g3d.θ[i,j,:] = g3d.θ[i,j,1] .+ (θ_1 .* g3d.σ) .+ (θ_2 .* g3d.σ.^2)
            else
                # set temperatures and grounded melt rate for cells without ice
                g3d.θ[i,j,:] .= params.default_temperature
                gh.θ_ave[i,j] = params.default_temperature_ave
                gh.grounded_basal_melt[i,j] = 0.
            end
        end
    end

    if minimum(gh.grounded_basal_melt[gh.mask]) < -1.0e-10
        println("WARNING: minimum grounded melt rate is negative")
    end

    if maximum(gh.θ_ave[gh.mask] .- (thermo_dynamics.kelvin_conversion .- params.density_ice .* params.g .* thermo_dynamics.melting_point_coeff .* gh.h[gh.mask]./2.)) > 0.001
        println("WARNING: maximum depth-averaged temperature is above pressure melting point")
    end

    return model
end