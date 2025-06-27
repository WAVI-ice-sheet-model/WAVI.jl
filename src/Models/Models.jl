module Models

using Parameters
using Setfield

using WAVI: AbstractField, AbstractGrid, AbstractMeltRate, AbstractModel, AbstractSpec
using WAVI.Fields
using WAVI.Grids
using WAVI.MeltRates
using WAVI.Parameters

export Model, update_state!

struct BasicSpec <: AbstractSpec 
    function BasicSpec()
        @info "Not implementing any parallel computations, running with BasicSpec"
        return new()
    end
end

struct Model{T,N,S,F,G,M} <: AbstractModel{T,N,S,F,G,M}
    grid    ::  G
    fields  ::  F
    params  ::  Params
    solver_params :: SolverParams
    spec   ::  S
    melt_rate :: M    
end

"""
    Model()

    Construct a WAVI.jl model object.

"""
function Model(grid::AbstractGrid, 
               bed_elevation::Union{Integer, Function}, 
               thickness::Integer,
               spec::AbstractSpec;
               initial_conditions::InitialConditions = InitialConditions(),
               params::Params = Params(),
               solver_params::SolverParams = SolverParams(),
               melt_rate = UniformMeltRate())

    bed_array = get_bed_elevation(bed_elevation, grid)
    
    # TODO: the passthrough of arguments like this is smelly - Configuration should be a type
    fields = GridField(grid, bed_array, thickness; initial_conditions, params, solver_params)
    model = Model{Real, Integer, BasicSpec, GridField, AbstractGrid, AbstractMeltRate}(grid, fields, Params(), SolverParams(), spec, melt_rate)
end

Model(grid, bed_elev, thickness; kw...) = Model(grid, bed_elev, thickness, BasicSpec(); kw...)
# This is to enable use of Setfield, which derives a parameter setup from the fields of an existing structure via JuliaObjects
Model(g, f, p, sp, s, m) = Model{Real, Integer, AbstractSpec, AbstractField, AbstractGrid, AbstractMeltRate}(g, f, p, sp, s, m)

include("utils.jl")

function Base.show(io::IO, m::Model)
    return print(io, "Characteristics of ", summary(m))
end

end
