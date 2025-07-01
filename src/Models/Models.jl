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

    Model(g, f, p, sp, s, m) = new{Float64,Int64,AbstractSpec,AbstractField,AbstractGrid,AbstractMeltRate}(g, f, p, sp, s, m)
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
               melt_rate = UniformMeltRate()) where {G<:AbstractGrid, S<:AbstractSpec}

    # FIXME: this all smells, hacking for threading
    bed_array = typeof(bed_elevation) <: AbstractArray ? bed_elevation : get_bed_elevation(bed_elevation, grid)
    
    # TODO: the passthrough of arguments like this is smelly - Configuration should be a type
    fields = GridField(grid, bed_array; initial_conditions, params, solver_params)
    model = Model(grid, fields, params, solver_params, spec, melt_rate)
    return model
end

Model(grid, bed_elev; kw...) = Model(grid, bed_elev, BasicSpec(); kw...)
Model(; grid, bed_elevation, spec, kw...) = Model(grid, bed_elevation, spec; kw...)

# This is to enable use of Setfield, which derives a parameter setup from the fields of an existing structure via JuliaObjects
# FIXME: this wasn't required in the original WAVI codebase. There is a constructor mapping issue somewhere.
#  If we specify T as Real, the solver will stop working - there is too much passthrough
#Model(g, f, p, sp, s, m) = Model{Float64,Int64}(g, f, p, sp, s, m)

include("utils.jl")

function Base.show(io::IO, m::Model)
    return print(io, "Characteristics of ", summary(m))
end

end
