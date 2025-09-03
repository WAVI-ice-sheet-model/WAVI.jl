module Fields

using LinearAlgebra
using LinearMaps
using Parameters
using Setfield          # TODO: InitialConditions using this, bit of an anti-pattern?
using SparseArrays

using WAVI: AbstractField, AbstractGrid
using WAVI.Grids
using WAVI.KroneckerProducts
using WAVI.Parameters
using WAVI.Utilities
using WAVI.Wavelets

export GridField, InitialConditions

include("UGrid.jl")
include("VGrid.jl")
include("HGrid.jl")
include("CGrid.jl")
include("SigmaGrid.jl")
include("InitialConditions.jl")
include("utils.jl")

"""
    Structure to hold all field variables in WAVI.jl
"""
struct GridField{T <: Real, N <: Integer} <: AbstractField{T, N}
    gh  :: HGrid{T,N}
    gu  :: UGrid{T,N}
    gv  :: VGrid{T,N}
    gc  :: CGrid{T,N}
    g3d :: SigmaGrid{T}

    wu  :: UWavelets{T,N}
    wv  :: VWavelets{T,N}    
end

"""
    setup_fields(grid, initial_conditions, solver_params, params, bed_array)

Acts as a constructor for the fields (no explicit constructor as fields `only ever called when setting up a model)
"""

function GridField(grid::AbstractGrid, bed_array;
                   initial_conditions::InitialConditions = InitialConditions(),
                   params::Params = Params(),
                   solver_params::SolverParams = SolverParams())

    initial_conditions = check_initial_conditions(initial_conditions, params, grid)

    ## Parameter fields checks 
    #if weertman c passed as a scalar, replace weertman_c parameters with matrix of this value
    if isa(params.weertman_c, Number) 
        params = @set params.weertman_c = params.weertman_c*ones(grid.nx,grid.ny)
    end
    #check size compatibility of resulting weertman C
    (size(params.weertman_c)==(grid.nx,grid.ny)) || throw(DimensionMismatch("Size of input weertman c ($(size(params.weertman_c))) must match grid size (i.e. $(grid.nx) x $(grid.ny))"))
    
    #if accumulation is passed as a scalar, replace accumulation parameters with matrix of this value
    if isa(params.accumulation_rate, Number) 
        params = @set params.accumulation_rate = params.accumulation_rate*ones(grid.nx,grid.ny)
    end
    #check size compatibility of resulting accumulation rate
    (size(params.accumulation_rate)==(grid.nx,grid.ny)) || throw(DimensionMismatch("Size of input accumulation ($(size(params.accumulation_rate))) must match grid size (i.e. $(grid.nx) x $(grid.ny))"))

    #if accumulation is passed as a scalar, replace accumulation parameters with matrix of this value
    if isa(params.glen_a_ref, Number) 
        params = @set params.glen_a_ref = params.glen_a_ref*ones(grid.nx,grid.ny)
    end
    #check size compatibility of resulting glen a ref
    (size(params.glen_a_ref)==(grid.nx,grid.ny)) || throw(DimensionMismatch("Size of input glen_a_ref ($(size(params.glen_a_ref))) must match grid size (i.e. $(grid.nx) x $(grid.ny))"))

    # TODO: grids are heavily reliant on the use of keyword arguments which do not allow specializations / multiple dispatch to work effectively

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
    grounded_fraction =  deepcopy(initial_conditions.initial_grounded_fraction)
    ηav = deepcopy(initial_conditions.initial_viscosity[:,:,1]) #set to the viscosity on the first level for now
    gh=HGrid(
        nxh=grid.nx,
        nyh=grid.ny,
        mask=h_mask,
        h_isfixed = grid.h_isfixed,
        b = bed_array,
        h = h,
        ηav = ηav,
        grounded_fraction = grounded_fraction
    )

    #u-grid
    gu=UGrid(
        nxu=grid.nx+1,
        nyu=grid.ny,
        dx=grid.dx,
        dy=grid.dy,
        mask=u_mask,
        u_isfixed=grid.u_isfixed,
        u=deepcopy(initial_conditions.initial_u_veloc),
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
        v=deepcopy(initial_conditions.initial_v_veloc),
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
        
    g3d=SigmaGrid(
        nxs=grid.nx,
        nys=grid.ny,
        nσs=grid.nσ,
        σ =grid.σ,
        η = η,
        θ = θ,
        Φ = Φ,
        glen_b = g3_glen_b,
        quadrature_weights = grid.quadrature_weights
    )

    #Wavelet-grid, u-component.
    wu=UWavelets(nxuw=grid.nx+1,nyuw=grid.ny,levels=solver_params.levels)

    #Wavelet-grid, v-component.
    wv=VWavelets(nxvw=grid.nx,nyvw=grid.ny+1,levels=solver_params.levels)
    return GridField(gh,gu,gv,gc,g3d,wu,wv)
end

end