"""
    update_tensile_strain_history!(model::AbstractModel)

Find the tensile strain energy history function from 'ice shelf' parts of strain rate tensor, neglecting all vertical shear.
"""
function update_tensile_strain_history!(model::AbstractModel{T,N}) where {T,N}
    @unpack gh,gu,gv,gc,g3d = model.fields
    @unpack params,solver_params=model

    tol = solver_params.tol_tensile_eig_solve
    maxiter = solver_params.maxiter_tensile_eig_solve
    max_strain_error = zero(T)

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
                    sigma3 = -params.density_ice.*params.g.*gh.h[i,j]*g3d.ζ[k]
                    sigma2 = resistive_stress_eig2[i,j] .+ sigma3
                    sigma1 = resistive_stress_eig1[i,j] .+ sigma3
                    degradation = one(T) - g3d.Φ[i,j,k] 
                    tensile_energy, strain_error = tensile_elastic_strain_energy(sigma1,sigma2,sigma3;
                         lambda=params.elastic_lambda,
                         mu=params.elastic_mu,
                         degradation=degradation,
                         tol=tol,
                         maxiter=maxiter) 
                    max_strain_error = max(strain_error,max_strain_error)
                    g3d.tensile_strain_history[i,j,k] = max(g3d.tensile_strain_history[i,j,k] , positivePart(
                       tensile_energy .- params.critical_elastic_energy))
                    phase_field = one(T) - one(T)./
                      (one(T) + 2*params.phase_field_length*g3d.tensile_strain_history[i,j,k]/params.energy_release_rate)
                    phase_field=clamp(phase_field,zero(T),one(T))
                    g3d.Φ[i,j,k]=one(T) - degradation_function(phase_field;k_reg = params.degradation_regularisation)    
                end
            end
        end
    end
    println(max_strain_error)
    return model
end

"""
    degradation_function(d)

Degradation function used in phase field model.
"""
function degradation_function(d;k_reg=0)
  degradation = k_reg + (one(k_reg)-k_reg)*(one(d)-d)^2
  return degradation
end

"""
    tensile_ealstic_strain_energy(sigma1,sigma2,sigma3)

Find the tensile part of the elastic strain energy from the three principal stresses using a fixed-point iteration
"""
function tensile_elastic_strain_energy(sigma1::T,sigma2::T,sigma3::T;lambda,mu,degradation=1.0,tol=eps(),maxiter=100) where T

    p = -(sigma1 .+ sigma2 .+ sigma3)./3.0
    K = lambda .+ (2/3)*mu
    eig3StartGuess = undegrade(-p./K,degradation) 
    eig3 = eig3StartGuess
    f3 = degrade(eig3,degradation)
    strain_error = Inf
    iterations=0
    alpha = 1.0

    while (strain_error > tol) & (iterations < maxiter)
        iterations += 1
        f3last=f3
        eig1=undegrade(f3+(sigma1-sigma3)./(2*mu),degradation)
        eig2=undegrade(f3+(sigma2-sigma3)./(2*mu),degradation)
        eig3=undegrade(f3,degradation)
        trace=eig1+eig2+eig3
        f3new = (sigma3 - lambda.*degrade(trace,degradation))./(2*mu)
        degTrace = (trace >= 0) ? degradation : one(degradation)
        invDeg1 = (eig1 >= 0) ? one(degradation)./degradation : one(degradation)
        invDeg2 = (eig2 >= 0) ? one(degradation)./degradation : one(degradation)
        invDeg3 = (eig3 >= 0) ? one(degradation)./degradation : one(degradation)
        chi = (lambda./(2*mu)).*degTrace.*(invDeg1.+invDeg2.+invDeg3)
        r = alpha .* one(chi)./(one(chi) .+ chi)
        f3 = (one(r) .- r).*f3last .+ r.*f3new

        last_strain_error = strain_error 
        strain_error = maximum(abs.(f3last - f3new))
        if strain_error > last_strain_error
          alpha = alpha .* 0.5
          f3 = f3last
        end
    end


    eig1=undegrade(f3+(sigma1-sigma3)./(2*mu),degradation)
    eig2=undegrade(f3+(sigma2-sigma3)./(2*mu),degradation)
    eig3=undegrade(f3,degradation)
    trace = eig1+eig2+eig3

    positiveStrainEnergy=(lambda/2).*(positivePart(eig1+eig2+eig3)).^2+mu.*((positivePart(eig1)).^2+(positivePart(eig2)).^2+(positivePart(eig3)).^2)

    if strain_error > tol
      println(degradation," ",iterations," ",strain_error," ",positiveStrainEnergy)
    end

    return positiveStrainEnergy,strain_error
end

positivePart(x) = 0.5*(x.+abs.(x))
negativePart(x) = 0.5*(x.-abs.(x))
degrade(x,degradation) = degradation.*positivePart(x) .+ negativePart(x)
undegrade(x,degradation) = positivePart(x)./degradation .+ negativePart(x)