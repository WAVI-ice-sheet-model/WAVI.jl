export BasicParallelSpec, update_preconditioner!, precondition!

struct BasicParallelSpec <: AbstractParallelSpec end

"""
update_preconditioner!(model::AbstractModel,::BasicParallelSpec)

Update the preconditioner. For Basic Parallel Specification no action is needed. 

"""
function update_preconditioner!(model::AbstractModel,::BasicParallelSpec)
    return model
end

function precondition!(model::AbstractModel, ::BasicParallelSpec)
    @unpack solver_params=model

    x = WAVI.get_start_guess(model)  
    op = WAVI.get_op(model)
    b = WAVI.get_rhs(model)
    resid = WAVI.get_resid(x,op,b)
    WAVI.set_residual!(model,resid)
    rel_resid = norm(resid)/norm(b)
    converged = rel_resid < solver_params.tol_picard

    correction = zero(x)

    if ! converged
      p=WAVI.get_preconditioner(model,op)
      WAVI.apply_precondition!(correction, p, resid)
      correction_coarse = WAVI.get_correction_coarse(p)
      WAVI.set_correction_coarse!(model,correction_coarse)
    end
    x .= x .+ correction
    WAVI.set_velocities!(model,x)

    return converged, rel_resid
end