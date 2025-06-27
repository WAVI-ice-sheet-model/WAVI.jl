export Preconditioner

#Struct to hold information about wavelet-based multigrid preconditioner.
@with_kw struct Preconditioner{T <: Real, 
                               O <:MapOrMatrix{T}, 
                               C <: MapOrMatrix{T}, 
                               R <: MapOrMatrix{T}, 
                               P <: MapOrMatrix{T}} <: AbstractPreconditioner{T,O,C,R,P}
    op::O
    op_diag::Vector{T} = diag(sparse(op))
    nsmooth::Integer = 5
    sweep::Vector{Integer}
    sweep_order::Vector{Integer} = unique(sweep)
    smoother_omega::T = 1.0
    restrict::R
    prolong::P
    op_coarse::C = restrict*op*prolong
    correction_coarse::Vector{T} = zeros(T,size(op_coarse,2))
    tol_coarse::T = 1e-7
    maxiter_coarse::Integer = 1000
end
