

struct ConstantDamage <: AbstractFracture end
struct DruckerPragerPhaseField <: AbstractFracture end

get_fracture(model::AbstractModel{T,N}) where {T,N} = model.fracture

"""
    update_damage!(model::AbstractModel)

Update the damage field and the strain history from 'ice shelf' parts of strain rate tensor, neglecting all vertical shear.
"""
update_damage!(model::AbstractModel{T,N};kwargs ...) where {T,N} = update_damage!(get_fracture(model),model; kwargs ...)

function update_damage!(fracture::ConstantDamage,model::AbstractModel{T,N};kwargs...) where {T,N}
  return model
end

function update_damage!(fracture::DruckerPragerPhaseField,model::AbstractModel{T,N}; store_strain_history = false) where {T,N}
    @unpack gh,gu,gv,gc,g3d = model.fields
    @unpack params,solver_params=model

    A_drucker_prager = 2.0* (params.ice_tensile_strength .* params.ice_compressive_strength) ./
                   (sqrt(3.0)*(params.ice_tensile_strength .+ params.ice_compressive_strength))
    B_drucker_prager = (params.ice_tensile_strength .- params.ice_compressive_strength) ./
                   (sqrt(3.0)*(params.ice_tensile_strength .+ params.ice_compressive_strength))

    resistive_stress_eig1=zeros(T,gh.nxh,gh.nyh)
    resistive_stress_eig2=zeros(T,gh.nxh,gh.nyh)

    resistive_stress_eig1[:] .= 
      gh.crop*gh.ηav[:].*( 3*(gu.∂x*(gu.crop*gu.u[:]) .+ gv.∂y*(gv.crop*gv.v[:])) .+
        sqrt.((gu.∂x*(gu.crop*gu.u[:]) .+ gv.∂y*(gv.crop*gv.v[:])).^2 .+
              (gc.cent*(gc.crop*( gu.∂y*(gu.crop*gu.u[:]) .+ gv.∂x*(gv.crop*gv.v[:])) )).^2
                 ) )

    resistive_stress_eig2[:] .= 
      gh.crop*gh.ηav[:].*( 3*(gu.∂x*(gu.crop*gu.u[:]) .+ gv.∂y*(gv.crop*gv.v[:])) .-
        sqrt.((gu.∂x*(gu.crop*gu.u[:]) .+ gv.∂y*(gv.crop*gv.v[:])).^2 .+
              (gc.cent*(gc.crop*( gu.∂y*(gu.crop*gu.u[:]) .+ gv.∂x*(gv.crop*gv.v[:])) )).^2
                 ) )

    for k=1:g3d.nσs
        for j=1:g3d.nys
            for i=1:g3d.nxs
                if gh.mask[i,j]
                    water_depth = max(params.sea_level_wrt_geoid - gh.s[i,j] + gh.h[i,j]*g3d.ζ[k],zero(T))
                    water_pressure = params.density_ocean.*params.g * water_depth
                    sigma3_effective = -params.density_ice.*params.g.*gh.h[i,j]*g3d.ζ[k] + water_pressure
                    sigma2_effective = resistive_stress_eig2[i,j] .+ sigma3_effective
                    sigma1_effective = resistive_stress_eig1[i,j] .+ sigma3_effective
                    degradation = one(T) - g3d.Φ[i,j,k] 
                    drucker_prager_elastic_energy = drucker_prager_elastic_strain_energy(sigma1_effective,sigma2_effective,sigma3_effective;
                         lambda=params.elastic_lambda,
                         mu=params.elastic_mu,
                         A = A_drucker_prager,
                         B = B_drucker_prager,
                         degradation = degradation) 
                    strain_history = max(g3d.strain_history[i,j,k] , drucker_prager_elastic_energy)
                    phase_field = one(T) - one(T)./
                      (one(T) + 2*params.phase_field_length*strain_history/params.energy_release_rate)
                    phase_field=clamp(phase_field,zero(T),one(T))
                    g3d.Φ[i,j,k]=one(T) - degradation_function(phase_field;k_reg = params.degradation_regularisation)  
                    store_strain_history && (g3d.strain_history[i,j,k] = strain_history)
                end
            end
        end
    end
    return model
end

"""
    update_strain_history!(model::AbstractModel)

Advect the strain history neglecting all vertical shear.
"""
function update_strain_history!(model::AbstractModel{T,N}) where {T,N} 
@unpack g3d = model.fields

update_damage!(model;store_strain_history = true)
advect3D!(g3d.strain_history,model)

return model

end

"""
    degradation_function(d)

Degradation function used in phase field model.
"""
function degradation_function(d::T;k_reg=zero(T)) where T
  degradation = k_reg + (one(k_reg)-k_reg)*(one(d)-d)^2
  return degradation
end

"""
    drucker_prager_elastic_strain_energy(sigma1,sigma2,sigma3)

Find the degraded part of the elastic strain energy from the three principal stresses.
"""
function drucker_prager_elastic_strain_energy(sigma1::T,sigma2::T,sigma3::T;lambda,mu,A,B,degradation=one(T)) where T

    p = -(sigma1 .+ sigma2 .+ sigma3)./3.0
    K = lambda .+ (2/3)*mu

    tau_deviatoric = sqrt.( 0.5*( (sigma1 .+ p).^2 .+ (sigma2 .+ p).^2 .+ (sigma3 .+ p).^2 ) ) 

    p_offset = p .- A./(3.0 .* B)

    if tau_deviatoric >= -3.0 * B * p_offset
       if tau_deviatoric >= mu*p_offset/(3.0*B*K)
          positiveStrainEnergy = (3.0*B*p_offset + tau_deviatoric).^2 ./ ( (2.0 *mu .+ 18.0 * B.^2 .* K) .* degradation.^2)
       else
          positiveStrainEnergy = ( p_offset.^2 ./ (2.0 * K) .+ tau_deviatoric.^2 ./ (2.0 *mu) ) ./ (degradation.^2) 
       end
    else
      positiveStrainEnergy = zero(T)
    end
    return positiveStrainEnergy
end

positivePart(x) = 0.5*(x.+abs.(x))
negativePart(x) = 0.5*(x.-abs.(x))
degrade(x,degradation) = degradation.*positivePart(x) .+ negativePart(x)
undegrade(x,degradation) = positivePart(x)./degradation .+ negativePart(x)