module Models

using Parameters
using Setfield

using WAVI: AbstractField, AbstractGrid, AbstractMeltRate, AbstractModel, AbstractSpec
using WAVI.Fields
using WAVI.Grids
using WAVI.MeltRates
using WAVI.Parameters

export Model, update_state!

"""
Struct to represent the basic specification for a model

"""
struct BasicSpec <: AbstractSpec 
    function BasicSpec()
        @info "Not implementing any parallel computations, running with BasicSpec"
        return new()
    end
end

struct Model{T,N,S,F,G,M<:AbstractMeltRate} <: AbstractModel{T,N,S,F,G}
    grid    ::  G
    fields  ::  F
    params  ::  Params
    solver_params :: SolverParams
    spec   ::  S
    melt_rate :: M

    Model{T,N,S,F,G,M}(g, f, p, sp, s, m) where {T,N,S,F,G,M} = new{T,N,S,F,G,M}(g, f, p, sp, s, m)
end

"""
    Model()

    Construct a WAVI.jl model object.

"""
function Model(grid::G, 
               bed_elevation::Union{Integer, Function, AbstractArray}, 
               spec::S;
               initial_conditions::InitialConditions = InitialConditions(),
               params::Params = Params(),
               solver_params::SolverParams = SolverParams(),
               melt_rate::M = UniformMeltRate()) where {G<:AbstractGrid, S<:AbstractSpec, M<:AbstractMeltRate}

    # FIXME: this all smells, hacking for threading
    bed_array = typeof(bed_elevation) <: AbstractArray ? bed_elevation : get_bed_elevation(bed_elevation, grid)
    
    # TODO: the passthrough of arguments like this is smelly - Configuration should be a type
    fields = GridField(grid, bed_array; initial_conditions, params, solver_params)
    model = Model{Float64, Int64, S, GridField, G, M}(grid, fields, params, solver_params, spec, melt_rate)
    return model
end

Model(grid, bed_elev; kw...) = Model(grid, bed_elev, BasicSpec(); kw...)
Model(; grid, bed_elevation, spec, kw...) = Model(grid, bed_elevation, spec; kw...)

# This is to enable use of Setfield, which derives a parameter setup from the fields of an existing structure via JuliaObjects
# FIXME: this wasn't required in the original WAVI codebase. 
#   Model needs to be in a position to have type analysis on it's properties to recreate the instance of it by ConstructionBase, which requires some refactoring
#   Ref: https://juliaobjects.github.io/ConstructionBase.jl/dev/#type-tips
Model(g::G, f::F, p::P, sp::SP, s::S, m::M) where {G<:AbstractGrid,F<:AbstractField,P<:Params,SP<:SolverParams,S<:AbstractSpec,M<:AbstractMeltRate} = 
    Model{Float64,Int64,S,F,G,M}(g, f, p, sp, s, m)

##
# Global domain alterations
#
Base.propertynames(model::Model{T,N,S,F,G,M}, private::Bool) where {T,N,S,F,G,M} = (fieldnames(typeof(model)..., :global_fields))

# FIXME: this is a sign of a frustration in WAVIs structural layout - too many deep nested structures accessed through high level passing
#  which inhibits multiple dispatch
function Base.getproperty(model::Model{T,N,S,F,G,M}, s::Symbol) where {T,N,S,F,G,M}
    if s == :global_fields
        return getfield(model, :fields)
    elseif s == :global_grid
        return getfield(model, :grid)
    end
    getfield(model, s)
end


include("utils.jl")

function Base.show(io::IO, m::Model)
    return print(io, "Characteristics of ", summary(m))
end

end
