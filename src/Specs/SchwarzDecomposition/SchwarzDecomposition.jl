using WAVI
import WAVI: AbstractModel
using WAVI.Models: Model

using Parameters
using Setfield

export  schwarzModel, 
        schwarzRestrictVelocities!, 
        schwarzProlongVelocities!,
        schwarzPartitionOfUnity




"""
schwarzModel(model::AbstractModel;igrid=1,jgrid=1,ngridsx=1,ngridsy=1,overlap=1)

Creates a subdomain model (model_g) used for parallel Schwarz iterations.
The necessary quantities are transferred from the full domain model.

Inputs:
model:        Full domain model
igrid:        Cartesian index of subdomain in x-direction.
jgrid:        Cartesian index of subdomain in y-direction.
ngridsx:      Number of subdomains in x-direction.
ngridsy:      Number of subdomains in y-direction.
overlap:      Overlap of subdomains on h-grid.

Output:
model_g:      Subdomain model

"""
function schwarzModel(model::AbstractModel;igrid=1,jgrid=1,ngridsx=1,ngridsy=1,overlap=1)
    @unpack nx,ny,dx,dy,nσ,x0,y0,h_mask,h_isfixed,hyd_potential_isfixed,u_iszero,v_iszero,u_isfixed,v_isfixed,quadrature_weights,σ = model.grid
    @unpack gh,gu,gv,g3d = model.fields
    
    @assert rem(nx,ngridsx)==0 "Model domain is not an integer number of subdomains in x-direction"
    @assert rem(ny,ngridsy)==0 "Model domain is not an integer number of subdomains in y-direction"
    (1 <= igrid <= ngridsx) || throw(BoundsError(model,igrid))
    (1 <= jgrid <= ngridsy) || throw(BoundsError(model,jgrid))

    nx_domain=div(nx,ngridsx)
    ny_domain=div(ny,ngridsy)
    i_start_domain=(igrid-1)*nx_domain+1
    i_stop_domain=igrid*nx_domain
    j_start_domain=(jgrid-1)*ny_domain+1
    j_stop_domain=jgrid*ny_domain

    i_start_g = max(i_start_domain - overlap, 1)
    i_stop_g = min(i_stop_domain + overlap, nx)
    j_start_g = max(j_start_domain - overlap, 1)
    j_stop_g = min(j_stop_domain + overlap, ny)

    nx_g = i_stop_g - i_start_g + 1
    ny_g = j_stop_g - j_start_g + 1
    dx_g = dx
    dy_g = dy
    nσ_g = nσ
    x0_g = x0 + (i_start_g-1)*dx
    y0_g = y0 + (j_start_g-1)*dy
    h_mask_g = h_mask[i_start_g:i_stop_g,j_start_g:j_stop_g]
    h_isfixed_g = h_isfixed[i_start_g:i_stop_g,j_start_g:j_stop_g]
    hyd_potential_isfixed_g = hyd_potential_isfixed[i_start_g:i_stop_g,j_start_g:j_stop_g]

    u_iszero_g = u_iszero[i_start_g:i_stop_g+1,j_start_g:j_stop_g]
    v_iszero_g = v_iszero[i_start_g:i_stop_g,j_start_g:j_stop_g+1]

    u_isfixed_g = u_isfixed[i_start_g:i_stop_g+1,j_start_g:j_stop_g] 
    v_isfixed_g = v_isfixed[i_start_g:i_stop_g,j_start_g:j_stop_g+1]

    #Set halo points as fixed velocity points
    (igrid==1) || (u_isfixed_g[1,:] .= true)
    (igrid==ngridsx) || (u_isfixed_g[nx_g+1,:] .= true)
    (igrid==1) || (v_isfixed_g[1,:] .= true)
    (igrid==ngridsx) || (v_isfixed_g[nx_g,:] .= true)
    (jgrid==1) || (v_isfixed_g[:,1] .= true)
    (jgrid==ngridsy) || (v_isfixed_g[:,ny_g+1] .= true)
    (jgrid==1) || (u_isfixed_g[:,1] .= true)
    (jgrid==ngridsy) || (u_isfixed_g[:,ny_g] .= true)

    quadrature_weights_g = quadrature_weights 
    σ_g = σ 

    grid_g=Grid(
    nx=nx_g,
    ny=ny_g,
    dx=dx_g,
    dy=dy_g,
    nσ=nσ_g,
    x0=x0_g,
    y0=y0_g,
    h_mask = h_mask_g,
    h_isfixed = h_isfixed_g,
    hyd_potential_isfixed = hyd_potential_isfixed_g,
    u_iszero = u_iszero_g,
    v_iszero = v_iszero_g,
    u_isfixed = u_isfixed_g,
    v_isfixed = v_isfixed_g,
    quadrature_weights = quadrature_weights_g,
    σ = σ_g)

    bed_elevation_g = model.fields.gh.b[i_start_g:i_stop_g,j_start_g:j_stop_g]

    params_g = Params(model.params,model.grid,(i_start_g,i_stop_g,j_start_g,j_stop_g))

    sliding_law_g = model.sliding_law
    if isdefined(sliding_law_g, :drag_coefficient)
        sliding_law_g = @set sliding_law_g.drag_coefficient = sliding_law_g.drag_coefficient[i_start_g:i_stop_g,j_start_g:j_stop_g]
    end

    basal_hydrology_g = model.basal_hydrology
    thermo_dynamics_g = model.thermo_dynamics

    solver_params_g=model.solver_params

    initial_thickness_g = gh.h[i_start_g:i_stop_g,j_start_g:j_stop_g]
    initial_grounded_fraction_g = gh.grounded_fraction[i_start_g:i_stop_g,j_start_g:j_stop_g]
    initial_u_veloc_g = gu.u[i_start_g:i_stop_g+1,j_start_g:j_stop_g]
    initial_v_veloc_g = gv.v[i_start_g:i_stop_g,j_start_g:j_stop_g+1]
    initial_viscosity_g = g3d.η[i_start_g:i_stop_g,j_start_g:j_stop_g,:]
    initial_temperature_g = g3d.θ[i_start_g:i_stop_g,j_start_g:j_stop_g,:]
    initial_damage_g = g3d.Φ[i_start_g:i_stop_g,j_start_g:j_stop_g,:]
    initial_strain_history_g = g3d.strain_history[i_start_g:i_stop_g,j_start_g:j_stop_g,:]
    initial_basal_water_thickness_g = gh.basal_water_thickness[i_start_g:i_stop_g,j_start_g:j_stop_g]
    initial_hydraulic_potential_b_g = gh.hydraulic_potential_b[i_start_g:i_stop_g,j_start_g:j_stop_g]
    initial_effective_pressure_g = gh.effective_pressure[i_start_g:i_stop_g,j_start_g:j_stop_g]
    initial_basal_melt_g = gh.basal_melt[i_start_g:i_stop_g,j_start_g:j_stop_g]
    initial_θ_ave_g = gh.θ_ave[i_start_g:i_stop_g,j_start_g:j_stop_g]
    initial_preBfactor_g = gh.preBfactor[i_start_g:i_stop_g,j_start_g:j_stop_g]

    initial_conditions_g=InitialConditions(
        initial_thickness = initial_thickness_g,
        initial_grounded_fraction = initial_grounded_fraction_g,
        initial_u_veloc = initial_u_veloc_g,
        initial_v_veloc = initial_v_veloc_g,
        initial_viscosity = initial_viscosity_g,
        initial_temperature = initial_temperature_g,
        initial_damage = initial_damage_g,
        initial_strain_history = initial_strain_history_g,
        initial_basal_water_thickness = initial_basal_water_thickness_g,
        initial_hydraulic_potential_b = initial_hydraulic_potential_b_g,
        initial_effective_pressure = initial_effective_pressure_g,
        initial_basal_melt = initial_basal_melt_g,
        initial_θ_ave = initial_θ_ave_g,
        initial_preBfactor = initial_preBfactor_g)

    shelf_melt_rate_g=model.shelf_melt_rate

    spec_g = BasicSpec()

    model_g=Model(
        grid = grid_g, 
        bed_elevation = bed_elevation_g,
        params = params_g,
        solver_params = solver_params_g,
        initial_conditions = initial_conditions_g,
        shelf_melt_rate = shelf_melt_rate_g,
        spec = spec_g,
        basal_hydrology = basal_hydrology_g,
        sliding_law = sliding_law_g,
        thermo_dynamics = thermo_dynamics_g)

    return model_g
end


"""
schwarzRestrictVelocities!(model_g::AbstractModel,model::AbstractModel;igrid=1,jgrid=1,ngridsx=1,ngridsy=1,overlap=1)

Schwarz restriction. Transfers velocities from the full domain (model) to the subdomain (model_g).
model_g:      Subdomain model
model:        Full domain model
igrid:        Cartesian index of subdomain in x-direction.
jgrid:        Cartesian index of subdomain in y-direction.
ngridsx:      Number of subdomains in x-direction.
ngridsy:      Number of subdomains in y-direction.
overlap:      Overlap of subdomains on h-grid.

"""
function schwarzRestrictVelocities!(model_g::AbstractModel, 
                                    model::AbstractModel; 
                                    igrid=1, jgrid=1, ngridsx=1, ngridsy=1, overlap=1)
    @unpack nx,ny = model.grid
    
    @assert rem(nx,ngridsx)==0 "Model domain is not an integer number of subdomains in x-direction"
    @assert rem(ny,ngridsy)==0 "Model domain is not an integer number of subdomains in y-direction"
    (1 <= igrid <= ngridsx) || throw(BoundsError(model,igrid))
    (1 <= jgrid <= ngridsy) || throw(BoundsError(model,jgrid))

    nx_domain=div(nx,ngridsx)
    ny_domain=div(ny,ngridsy)
    i_start_domain=(igrid-1)*nx_domain+1
    i_stop_domain=igrid*nx_domain
    j_start_domain=(jgrid-1)*ny_domain+1
    j_stop_domain=jgrid*ny_domain

    i_start_g = max(i_start_domain - overlap, 1)
    i_stop_g = min(i_stop_domain + overlap, nx)
    j_start_g = max(j_start_domain - overlap, 1)
    j_stop_g = min(j_stop_domain + overlap, ny)

    model_g.fields.gu.u .= model.fields.gu.u[i_start_g:i_stop_g+1,j_start_g:j_stop_g]
    model_g.fields.gv.v .= model.fields.gv.v[i_start_g:i_stop_g,j_start_g:j_stop_g+1]

    return model_g
end

"""
schwarzProlongVelocities!(model::AbstractModel,model_g::AbstractModel;igrid=1,jgrid=1,ngridsx=1,ngridsy=1,overlap=1,damping=0.0)

Schwarz prolongation. Transfers velocities from the subdomain (model_g) to the full domain (model)
model:        Full domain model
model_g:      Subdomain model
igrid:        Cartesian index of subdomain in x-direction.
jgrid:        Cartesian index of subdomain in y-direction.
ngridsx:      Number of subdomains in x-direction.
ngridsy:      Number of subdomains in y-direction.
overlap:      Overlap of subdomains on h-grid.
damping:      Damping factor.

"""
function schwarzProlongVelocities!(model::AbstractModel,
                                   model_g::AbstractModel;
                                   igrid=1, jgrid=1, ngridsx=1, ngridsy=1, overlap=1, damping=0.0)
    @unpack nx,ny = model.grid
    
    @assert rem(nx,ngridsx)==0 "Model domain is not an integer number of subdomains in x-direction"
    @assert rem(ny,ngridsy)==0 "Model domain is not an integer number of subdomains in y-direction"
    (1 <= igrid <= ngridsx) || throw(BoundsError(model,igrid))
    (1 <= jgrid <= ngridsy) || throw(BoundsError(model,jgrid))

    nx_domain=div(nx,ngridsx)
    ny_domain=div(ny,ngridsy)
    i_start_domain=(igrid-1)*nx_domain+1
    i_stop_domain=igrid*nx_domain
    j_start_domain=(jgrid-1)*ny_domain+1
    j_stop_domain=jgrid*ny_domain

    i_start_g = max(i_start_domain - overlap, 1)
    i_stop_g = min(i_stop_domain + overlap, nx)
    j_start_g = max(j_start_domain - overlap, 1)
    j_stop_g = min(j_stop_domain + overlap, ny)

    nx_g = i_stop_g - i_start_g + 1
    ny_g = j_stop_g - j_start_g + 1

    upou = partition_of_unity(nx_g+1,ny_g,igrid==1,igrid==ngridsx,jgrid==1,jgrid==ngridsy,2*overlap,2*overlap-1)
    vpou = partition_of_unity(nx_g,ny_g+1,igrid==1,igrid==ngridsx,jgrid==1,jgrid==ngridsy,2*overlap-1,2*overlap)

    model.fields.gu.u[i_start_g:i_stop_g+1,j_start_g:j_stop_g] .+= (one(damping)-damping) .* upou .* model_g.fields.gu.u
    model.fields.gv.v[i_start_g:i_stop_g,j_start_g:j_stop_g+1] .+= (one(damping)-damping) .* vpou .* model_g.fields.gv.v

    return model
end
