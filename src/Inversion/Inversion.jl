include("InversionDataHGrid.jl")
include("InversionDataUGrid.jl")
include("InversionDataVGrid.jl")
include("DataFields.jl")


struct Inversion{T <: Real, N <: Integer} 
    grid::Grid{T,N}
    data_fields::DataFields{T,N}
end

"""
    Inversion(;
        grid = nothing
        )

Construct a WAVI.jl inversion object.

Keyword arguments
=================
- `grid`: (required) an instance of a `Grid` object, which defines the computational grid

"""
function Inversion(;
    grid = nothing,
    speed_u = nothing,
    speed_u_mask = nothing,
    speed_v = nothing,
    speed_v_mask = nothing,
    dhdt = nothing,
    accumulation_rate = nothing,
    dhdtacc_mask = nothing)

    #check that a grid and bed has been inputted
    ~(grid === nothing) || throw(ArgumentError("You must specify an input grid"))

    #Setup the fields 
    data_fields = setup_datafields(grid,speed_u,speed_u_mask,speed_v,speed_v_mask,dhdt,accumulation_rate,dhdtacc_mask)

    #Use type constructor to build initial state with no extra physics
    inversion=Inversion(grid,data_fields)

    return inversion
end

include("run_inversion.jl")
#include("InversionData.jl")
#include("InversionDataUGrid.jl")
#include("InversionDataVGrid.jl")
#include("InversionDataHGrid.jl")
