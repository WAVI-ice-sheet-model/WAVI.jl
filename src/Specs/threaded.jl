export ThreadedSpec

using Parameters

using WAVI: AbstractModel, AbstractSpec

include("SchwarzDecomposition/SchwarzDecomposition.jl")
using .SchwarzDecomposition

import WAVI.Processes: update_preconditioner!, precondition!

"""
Struct to represent the shared memory parallel specification of a model.

"""
@with_kw struct ThreadedSpec{T<:Real, N<:Integer} <: AbstractSpec
    ngridsx::N = 1 
    ngridsy::N = 1
    overlap::N = 1
    niterations::N = 0
    damping::T = 0.0
    schwarzModelArray::Array{AbstractModel,2} = Array{AbstractModel,2}(undef,ngridsx,ngridsy)
end

"""
update_preconditioner!(model::AbstractModel, ::ThreadedSpec)

Update the preconditioner. 
For shared memory parallelism, update the subdomain models used by restricted additive Schwarz preconditioner (RAS).

"""
function update_preconditioner!(model::AbstractModel, ::ThreadedSpec)
    @unpack ngridsx, ngridsy, overlap = model.parallel_spec

    @sync for igrid = 1:ngridsx
        for jgrid = 1:ngridsy
            Threads.@spawn begin
                model.parallel_spec.schwarzModelArray[igrid,jgrid] = schwarzModel(model;
                                                                                  igrid=igrid,
                                                                                  jgrid=jgrid,
                                                                                  ngridsx=ngridsx,
                                                                                  ngridsy=ngridsy, 
                                                                                  overlap=overlap)
            end
        end
    end
end


"""
precondition!(model::AbstractModel,::ThreadedSpec)

Apply restricted additive Schwarz preconditioner (RAS) using shared memory parallelism.

"""
function precondition!(model::AbstractModel, ::ThreadedSpec)
    @unpack ngridsx, ngridsy, overlap, niterations, schwarzModelArray, damping = model.parallel_spec
    @unpack solver_params = model

    op = WAVI.get_op(model)
    b = WAVI.get_rhs(model)
    resid = WAVI.get_resid(x,op,b)
    WAVI.set_residual!(model,resid)
    rel_resid = norm(resid)/norm(b)
    converged = rel_resid < solver_params.tol_picard

    if ! converged
        for iteration = 1:niterations
            println("Schwarz iteration $iteration")
            @sync for igrid = 1:ngridsx
                for jgrid = 1:ngridsy
                    Threads.@spawn begin                
                        model_g = schwarzModelArray[igrid,jgrid]

                        schwarzRestrictVelocities!(
                            model_g::AbstractModel,
                            model::AbstractModel;
                            igrid=igrid,
                            jgrid=jgrid,
                            ngridsx=ngridsx,
                            ngridsy=ngridsy,
                            overlap=overlap)
                    end
                end
            end
            
            @sync for igrid = 1:ngridsx
                for jgrid = 1:ngridsy
                    Threads.@spawn begin
                        model_g = schwarzModelArray[igrid,jgrid]
                        WAVI.update_state!(model_g)
                    end
                end
            end

            model.fields.gu.u[:,:] .= damping .* model.fields.gu.u
            model.fields.gv.v[:,:] .= damping .* model.fields.gv.v
            
            threadLock=ReentrantLock()
            @sync for igrid = 1:ngridsx
                for jgrid = 1:ngridsy
                    Threads.@spawn begin
                        model_g = schwarzModelArray[igrid,jgrid]

                        lock(threadLock)
                        try
                            schwarzProlongVelocities!(model::AbstractModel,
                                model_g::AbstractModel;
                                igrid=igrid,
                                jgrid=jgrid,
                                ngridsx=ngridsx,
                                ngridsy=ngridsy,
                                overlap=overlap,
                                damping=damping)
                        finally
                            unlock(threadLock)
                        end
                    end                    
                end
            end
            println("")
        end

    end
    return converged, rel_resid
end
