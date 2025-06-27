"""
    update_tensile_strain_history!(model::AbstractModel)

Find the tensile strain energy history function from 'ice shelf' parts of strain rate tensor, neglecting all vertical shear.
"""
function update_tensile_strain_history!(model::AbstractModel{T,N}) where {T,N}
    @unpack gh,gu,gv,gc,g3d = model.fields
    @unpack params,solver_params=model

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

    advect3D!(g3d.tensile_strain_history,model)

    for k=1:g3d.nσs
        for j=1:g3d.nys
            for i=1:g3d.nxs
                if gh.mask[i,j]
                    viscosity_factor=g3d.η[i,j,k]/gh.ηav[i,j]
                    sigma3 = -params.density_ice.*params.g.*gh.h[i,j]*g3d.ζ[k]
                    sigma2 = viscosity_factor*resistive_stress_eig2[i,j] .+ sigma3
                    sigma1 = viscosity_factor*resistive_stress_eig1[i,j] .+ sigma3
                    degradation = one(T) - g3d.Φ[i,j,k] 
                    g3d.tensile_strain_history[i,j,k] = max(g3d.tensile_strain_history[i,j,k] , positivePart(
                        tensile_elastic_strain_energy(sigma1,sigma2,sigma3;lambda=params.elastic_lambda,mu=params.elastic_mu,degradation=degradation,tol=1.0e-7) .- params.critical_elastic_energy))
                    phase_field = one(T) - one(T)./
                      (one(T) + 2*params.phase_field_length*g3d.tensile_strain_history[i,j,k]/params.energy_release_rate)
                    g3d.Φ[i,j,k]=one(T) - degradation_function(phase_field;k_reg = params.degradation_regularisation)    
                end
            end
        end
    end

    return model
end

"""
    degradation_function(d)

Degradation function used in phase field model.
"""
degradation_function(d;k_reg=0) = k_reg + (one(k_reg)-k_reg)*(one(d)-d)^2


"""
    tensile_ealstic_strain_energy(sigma1,sigma2,sigma3)

Find the tensile part of the elastic strain energy from the three principal stresses using a fixed-point iteration
"""
function tensile_elastic_strain_energy(sigma1::T,sigma2::T,sigma3::T;lambda,mu,degradation=1.0,tol=eps()) where T

    eig3Start=zero(T)
    eig3=eig3Start
    f3 = degrade(eig3,degradation)
    err = Inf

    while err > tol
        eig1=undegrade(f3+(sigma1-sigma3)./(2*mu),degradation)
        eig2=undegrade(f3+(sigma2-sigma3)./(2*mu),degradation)
        trace=eig1+eig2+eig3
        f3last=f3
        f3 = 0.5*(f3+(sigma3 - lambda.*degrade(trace,degradation))./(2*mu))
        eig3=undegrade(f3,degradation)
        err = maximum(abs.(f3last - (sigma3 - lambda.*degrade(trace,degradation))./(2*mu)))
        println(err)
    end

    eig1=undegrade(f3+(sigma1-sigma3)./(2*mu),degradation)
    eig2=undegrade(f3+(sigma2-sigma3)./(2*mu),degradation)
    eig3=undegrade(f3,degradation)

    positiveStrainEnergy=(lambda/2).*(positivePart(eig1+eig2+eig3)).^2+mu.*((positivePart(eig1)).^2+(positivePart(eig2)).^2+(positivePart(eig3)).^2)

    return positiveStrainEnergy
end

positivePart(x) = 0.5*(x+abs(x))
negativePart(x) = 0.5*(x-abs(x))
degrade(x,degradation) = degradation.*positivePart(x) + negativePart(x)
undegrade(x,degradation) = positivePart(x)./degradation + negativePart(x)