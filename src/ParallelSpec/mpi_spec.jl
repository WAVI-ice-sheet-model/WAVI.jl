export MPISpec, update_preconditioner!, precondition!

using MPI

@with_kw struct MPISpec{T, N} <: AbstractParallelSpec
    ngridsx::N = 1 
    ngridsy::N = 1
    overlap::N = 1
    niterations::N = 0
    damping::T = 0.0
    mpiModelArray::Array{AbstractModel,2} = Array{AbstractModel,2}(undef,ngridsx,ngridsy)
end

update_preconditioner!(model::AbstractModel, ::MPISpec) = update_preconditioner!(model, BasicSpec())

precondition!(model::AbstractModel, ::MPISpec) = precondition!(model, BasicSpec())