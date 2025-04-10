export DistributedSpec, update_preconditioner!, precondition!

using .SchwarzDecomposition
using Distributed

"""
Struct to represent the shared memory parallel specification of a model.

"""
@with_kw struct DistributedSpec{T, N} <: AbstractParallelSpec
    ngridsx::N = 1 
    ngridsy::N = 1
    overlap::N = 1
    niterations::N = 0
    damping::T = 0.0
    distModelArray::Array{AbstractModel,2} = Array{AbstractModel,2}(undef,ngridsx,ngridsy)
    distProcessIds::Vector{Integer} = []
end

"""
update_preconditioner!(model::AbstractModel,::DistributedSpec)

Update the preconditioner. 

"""
function update_preconditioner!(model::AbstractModel, ::DistributedSpec)
    @unpack ngridsx, ngridsy, overlap, distProcessIds = model.parallel_spec
    
    # TODO: this won't work for non-LocalManager, revise.
    
    if length(distProcessIds) == 0
        ids = Distributed.addprocs(ngridsy; restrict=true)
        @info "Started processes $(ids)"
        append!(distProcessIds, ids)
    end

    # TODO: run loop across full x by y 
    for igrid = 1:ngridsx
        @sync @distributed for jgrid = 1:ngridsy
            model.parallel_spec.distModelArray[igrid,jgrid] = 
                schwarzModel(model;
                    igrid=igrid,
                    jgrid=jgrid,
                    ngridsx=ngridsx,
                    ngridsy=ngridsy, 
                    overlap=overlap)
        end
    end
end

"""
precondition!(model::AbstractModel,::SharedMemorySpec)

Apply restricted additive Schwarz preconditioner (RAS) using shared memory parallelism.

"""
function precondition!(model::AbstractModel, ::DistributedSpec)
    @unpack ngridsx, ngridsy, overlap, niterations, distModelArray, damping = model.parallel_spec
    @unpack solver_params = model

    x = WAVI.get_start_guess(model)  
    op = WAVI.get_op(model)
    b = WAVI.get_rhs(model)
    resid = WAVI.get_resid(x,op,b)
    WAVI.set_residual!(model,resid)
    rel_resid = norm(resid)/norm(b)
    converged = rel_resid < solver_params.tol_picard

    if ! converged
        for iteration = 1:niterations
            println("Schwarz iteration $iteration")
            for igrid = 1:ngridsx
                @sync @distributed for jgrid = 1:ngridsy
                    model_g = distModelArray[igrid,jgrid]

                    schwarzRestrictVelocities!(
                        model_g::AbstractModel,
                        model::AbstractModel;
                        igrid=igrid,
                        jgrid=jgrid,
                        ngridsx=ngridsx,
                        ngridsy=ngridsy,
                        overlap=overlap)
                
                    WAVI.update_state!(model_g)
                    model.fields.gu.u[:,:] .= damping .* model.fields.gu.u
                    model.fields.gv.v[:,:] .= damping .* model.fields.gv.v
                    
                    schwarzProlongVelocities!(
                        model::AbstractModel,
                        model_g::AbstractModel;
                        igrid=igrid,
                        jgrid=jgrid,
                        ngridsx=ngridsx,
                        ngridsy=ngridsy,
                        overlap=overlap,
                        damping=damping)
                end
            end
            println("")
        end
    end
    return converged, rel_resid
end
