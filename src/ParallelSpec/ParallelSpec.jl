module ParallelSpec

using LinearAlgebra, Parameters, Setfield
using WAVI


"""
get_parallel_spec(model::AbstractModel)

Function to return the parallel specification of a model.

"""
get_parallel_spec(model::AbstractModel) = model.parallel_spec

precondition!(model::AbstractModel) = precondition!(model,get_parallel_spec(model))

"""
update_preconditioner!(model)

Update the preconditioner for Basic Parallel Specification or Shared Memory Parallel Specification as appropriate.

"""
update_preconditioner!(model::AbstractModel) = update_preconditioner!(model::AbstractModel, get_parallel_spec(model::AbstractModel))

include("basic_spec.jl")
include("shared_memory_spec.jl")

end