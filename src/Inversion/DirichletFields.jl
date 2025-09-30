
"""
    Structure to hold all dirichlet field variables in WAVI.jl
"""
struct DirichletFields{T <: Real, N <: Real}
    gh::HGrid{T,N}
    gu::UGrid{T,N}
    gv::VGrid{T,N}
    gc::CGrid{T,N}
    g3d::SigmaGrid{T,N}   
end

"""
    setup_dirichlet_fields(grid)

Acts as a constructor for the dirichletfields (no explicit constructor as dirichletfields `only ever called when setting up a model)
"""

function setup_dirichletfields(grid,bed_array,solver_params,initial_conditions,params)
    #Define masks for points on h-, u-, v- and c-grids that lie in model domain.
    h_mask = grid.h_mask 
    u_mask = get_u_mask(h_mask)
    v_mask = get_v_mask(h_mask)
    c_mask = get_c_mask(h_mask)

    #Remove all points on u- and v-grids with homogenous Dirichlet conditions.
    u_mask[grid.u_iszero].=false
    v_mask[grid.v_iszero].=false

    #h-grid
    #gh=HGrid(grid, params) #uncomment if using the explicit constructor method
      h =  deepcopy(initial_conditions.initial_thickness)
    gh=HGrid(
    nxh=grid.nx,
    nyh=grid.ny,
    mask=h_mask,
    h_isfixed = grid.h_isfixed,
    b = bed_array,
    h = h
   # ηav = ηav,
    #grounded_fraction = grounded_fraction
    )

    #u-grid
    gu=UGrid(
    nxu=grid.nx+1,
    nyu=grid.ny,
    dx=grid.dx,
    dy=grid.dy,
    mask=u_mask,
    u_isfixed=grid.u_isfixed,
  #  u=deepcopy(initial_conditions.initial_u_veloc),
    levels=solver_params.levels
    )

    #v-grid
    gv=VGrid(
    nxv=grid.nx,
    nyv=grid.ny+1,
    dx=grid.dx,
    dy=grid.dy,
    mask=v_mask,
    v_isfixed=grid.v_isfixed,
  #  v=deepcopy(initial_conditions.initial_v_veloc),
    levels=solver_params.levels
    )

    #c-grid
    gc=CGrid(
    nxc=grid.nx-1,
    nyc=grid.ny-1,
    mask=c_mask
    )

    #3D-grid
    η = deepcopy(initial_conditions.initial_viscosity)
    θ = deepcopy(initial_conditions.initial_temperature)
    Φ = deepcopy(initial_conditions.initial_damage)
    g3_glen_b = zeros(size(η))
    for i = 1:grid.nx
        for j = 1:grid.ny
            for k = 1:grid.nσ
               g3_glen_b[i,j,k] = glen_b.(θ[i,j,k],Φ[i,j,k],params.glen_a_ref[i,j], params.glen_n, params.glen_a_activation_energy, params.glen_temperature_ref, params.gas_const)
            end
        end
    end
      
    preBfactor=ones(grid.nx,grid.ny,grid.nσ)

    g3d=SigmaGrid(
    nxs=grid.nx,
    nys=grid.ny,
    nσs=grid.nσ,
    σ =grid.σ,
    η = η,
    θ = θ,
    Φ = Φ,
    glen_b = g3_glen_b,
    quadrature_weights = grid.quadrature_weights,
    preBfactor= preBfactor
    ) 

    return DirichletFields(gh,gu,gv,gc,g3d)
end
