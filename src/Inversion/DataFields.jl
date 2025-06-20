#include("InversionDataHGrid.jl")
#include("InversionDataUGrid.jl")
#include("InversionDataVGrid.jl")


"""
    Structure to hold all datafield variables in WAVI.jl
"""
struct DataFields{T <: Real, N <: Real}
    ghdata::InversionDataHGrid{T,N}
    gudata::InversionDataUGrid{T,N}
    gvdata::InversionDataVGrid{T,N}
end

"""
    setup_datafields(grid)

Acts as a constructor for the datafields (no explicit constructor as datafields `only ever called when setting up a model)
"""

function setup_datafields(grid,speed_u,speed_u_mask,speed_v,speed_v_mask,dhdt,accumulation_rate, dhdtacc_mask, model)
    #Define masks for points on h-, u-, v- and c-grids that lie in model domain.
    h_mask = dhdtacc_mask 
    u_mask = speed_u_mask
    v_mask = speed_v_mask

    #Remove all points on u- and v-grids with homogenous Dirichlet conditions.
    u_mask[grid.u_iszero].=false
    v_mask[grid.v_iszero].=false

    #h-grid
    #gh=HGrid(grid, params) #uncomment if usexing the explicit constructor method
    ghdata=InversionDataHGrid(
    nxh=grid.nx,
    nyh=grid.ny,
    mask=h_mask,
    h_isfixed = grid.h_isfixed,
    dhdt = dhdt,
    accumulation_rate = accumulation_rate
    )

    #u-grid
    gudata=InversionDataUGrid(
    nxu=grid.nx+1,
    nyu=grid.ny,
    mask=u_mask,
    u_isfixed=grid.u_isfixed,
    speed_u=speed_u,
    model=model
    )

    #v-grid
    gvdata=InversionDataVGrid(
    nxv=grid.nx,
    nyv=grid.ny+1,
    mask=v_mask,
    v_isfixed=grid.v_isfixed,
    speed_v=speed_v,
    model=model
    )

    return DataFields(ghdata,gudata,gvdata)
end
