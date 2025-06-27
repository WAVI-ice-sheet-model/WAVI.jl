using LinearAlgebra
using Parameters

using WAVI: AbstractModel

export update_preconditioner!, precondition!

function update_preconditioner!(model::AbstractModel)

end

function precondition!(model::AbstractModel)
    @unpack solver_params=model

    x = get_start_guess(model)  
    op = get_op(model)
    b = get_rhs(model)
    resid = get_resid(x,op,b)
    set_residual!(model,resid)
    rel_resid = norm(resid)/norm(b)
    converged = rel_resid < solver_params.tol_picard

    # FIXME: this is to provide some variation
    correction = randn(size(x)) / 500
    # correction = zero(x)

    if ! converged
      # FIXME: we ditch out of wavelet based preconditioning for MiniWAVI
      # p=get_preconditioner(model, op)
      # apply_precondition!(correction, p, resid)
      # correction_coarse = get_correction_coarse(p)
      # set_correction_coarse!(model,correction_coarse)
    end
    x .= x .+ correction
    set_velocities!(model,x)

    return converged, rel_resid
end
