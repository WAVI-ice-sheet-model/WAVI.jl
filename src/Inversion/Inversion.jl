include("InversionDataHGrid.jl")
include("InversionDataUGrid.jl")
include("InversionDataVGrid.jl")
include("DataFields.jl")
include("DirichletFields.jl")

#struct Inversion{T <: Real, N <: Integer} 
struct Inversion{T <: Real, N <: Integer, A, W, G, M <:AbstractMeltRate, PS <: AbstractParallelSpec} <: AbstractModel{T,N,M,PS}
    grid::Grid{T,N}
    inversion_params::InversionParams{T,A,W,G}
    data_fields::DataFields{T,N}
    fields::DirichletFields{T,N}
    melt_rate::M
    parallel_spec::PS
    inversion_output::InversionOutput{T}
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
    bed_elevation = nothing,
    inversion_params = InversionParams(),
    solver_params = SolverParams(),
    speed_u = nothing,
    speed_u_mask = nothing,
    speed_v = nothing,
    speed_v_mask = nothing,
    dhdt = nothing,
    accumulation_rate = nothing,
    dhdtacc_mask = nothing,
    melt_rate = UniformMeltRate(),
    parallel_spec = BasicParallelSpec(),
    initial_conditions = InitialConditions(),
    params = Params(),
    inversion_output = InversionOutput())

    #check that a grid and bed has been inputted
    ~(grid === nothing) || throw(ArgumentError("You must specify an input grid"))

    #Setup the fields 
    data_fields = setup_datafields(grid,speed_u,speed_u_mask,speed_v,speed_v_mask,dhdt,accumulation_rate,dhdtacc_mask)


    bed_array = zeros(grid.nx,grid.nx) #initialize a bed array
    try
    bed_array = get_bed_elevation(bed_elevation, grid)
    catch
    throw(ArgumentError("bed elevation must be of type function or array"))
    end

     #if glen_a_ref is passed as a scalar, replace glen_a_ref parameters with matrix of this value
     if isa(params.glen_a_ref, Number) 
        params = @set params.glen_a_ref = params.glen_a_ref*ones(grid.nx,grid.ny)
    end
    #check size compatibility of resulting glen a ref
    (size(params.glen_a_ref)==(grid.nx,grid.ny)) || throw(DimensionMismatch("Size of input glen_a_ref must match grid size (i.e. $(grid.nx) x $(grid.ny))"))



    #Setup the fields 
    fields = setup_dirichletfields(grid,bed_array,solver_params,initial_conditions,params)

    #Use type constructor to build initial state with no extra physics
    inversion=Inversion(grid,inversion_params,data_fields,fields,melt_rate,parallel_spec,inversion_output)

    return inversion
end

#include("run_inversion.jl")
include("inversion_utilities.jl")
#include("InversionData.jl")
#include("InversionDataUGrid.jl")
#include("InversionDataVGrid.jl")
#include("InversionDataHGrid.jl")
