#include("DirichletFields.jl")

#= struct Inversion{T <: Real, N <: Integer, M <:AbstractMeltRate, PS <: AbstractParallelSpec} <: AbstractModel{T,N,M,PS}
   # grid::Grid{T,N}
    model:: AbstractModel{T,N,M,PS}
    inversion_params::InversionParams{T}
    data::DataFields{T,N}
 #   fields::DirichletFields{T,N}
  #  melt_rate::M
  #  parallel_spec::PS
    inversion_output::InversionOutput{T}
end =#

struct Inversion{T <: Real, N <: Integer, M <: AbstractMeltRate, PS <: AbstractParallelSpec}
    model::AbstractModel{T,N,M,PS}
    inversion_params::InversionParams{T}
    data::DataFields{T,N}
    inversion_output::InversionOutput{T}
end

"""
        Inversion(;
            model = Model()
            inversion_params = InversionParams(),
            data = data,
            inversion_output = InversionOutput(),
            )

Construct a WAVI.jl inversion object.

Keyword arguments
=================
    - `model`: model
    - `inversion_params`: inversion parameters
    - `data`: data to be used in the inversion
    - `inversion_output`: output from the inversion

"""
function Inversion(;
    model = Model(),
    inversion_params = InversionParams(),
    data = DataFields(),
    inversion_output = InversionOutput())
   
   # grid = deepcopy(model.grid)
   # params=deepcopy(model.params)
   # melt_rate=deepcopy(model.melt_rate)
   # parallel_spec=deepcopy(model.parallel_spec)
   # solver_params=deepcopy(model.solver_params)

    model = deepcopy(model)

    if !all(data.ghdata.accumulation_rate .== model.params.accumulation_rate)
    println("WARNING: a different accumulation_rate is being used in model and inversion.data!")
    else
    println("The same accumulation_rate is being used in model and inversion.data!")
    end

   # bed_array = deepcopy(model.fields.gh.b)

    #Check initial conditions, and revert to default values if not
   # initial_conditions = check_initial_conditions(model.initial_conditions, model.params, model.grid)

    #Setup the fields 
  #  fields = setup_dirichletfields(grid,bed_array,solver_params,initial_conditions,params)

    #Use type constructor to build initial state with no extra physics but containing observational data
  #  inversion=Inversion(grid,inversion_params,data,fields,melt_rate,parallel_spec, inversion_output)
    inversion=Inversion(model,inversion_params,data,inversion_output)

    return inversion
end

include("inversion_utilities.jl")
