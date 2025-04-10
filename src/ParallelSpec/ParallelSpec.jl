module ParallelSpec

using LinearAlgebra, Parameters, Setfield
using WAVI

include("SchwarzDecomposition/SchwarzDecomposition.jl")

"""
get_parallel_spec(model::AbstractModel)

Function to return the parallel specification of a model, used for multiple dispatch

"""
get_parallel_spec(model::AbstractModel) = model.parallel_spec

precondition!(model::AbstractModel) = precondition!(model, get_parallel_spec(model))
update_preconditioner!(model::AbstractModel) = update_preconditioner!(model::AbstractModel, get_parallel_spec(model::AbstractModel))

include("basic_spec.jl")
include("shared_memory_spec.jl")
include("distributed_spec.jl")
include("mpi_spec.jl")
include("slurm_spec.jl")

end