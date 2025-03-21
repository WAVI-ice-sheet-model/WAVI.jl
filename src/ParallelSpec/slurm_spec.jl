export SlurmSpec, update_preconditioner!, precondition!

@with_kw struct SlurmSpec{T, N} <: AbstractParallelSpec
    ngridsx::N = 1 
    ngridsy::N = 1
    overlap::N = 1
    niterations::N = 0
    damping::T = 0.0
    nodeArray::Array{AbstractModel,2} = Array{AbstractModel,2}(undef,ngridsx,ngridsy)
end

update_preconditioner!(model::AbstractModel, ::SlurmSpec) = update_preconditioner!(model, BasicSpec())

precondition!(model::AbstractModel, ::SlurmSpec) = precondition!(model, BasicSpec())