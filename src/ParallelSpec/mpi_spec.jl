export MPISpec, update_preconditioner!, precondition!, update_velocities!

"""
VERY naive implementation of a workflow for coordination of domain decomposition solving 
across MPI workers 

TODO: What we really need to do is initialise the global model split across members
TODO: this might (?) be best achieved leveraging uniform buffers
TODO: Currently, we're looking at sequential, per-iteration transfer of velocities
TODO: and a copy per process of the global model (which at present is probably out of date each time)

See GH#91 for high level tasks

"""

using MPI
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
comm_size = MPI.Comm_size(comm)
root_rank = 0

@with_kw struct MPISpec{T <: Real, N <: Integer} <: AbstractParallelSpec
    ngridsx::N = 1 
    ngridsy::N = 1
    overlap::N = 0
    niterations::N = 0
    damping::T = 0.0
    mpiModelArray::Array{AbstractModel,2} = Array{AbstractModel,2}(undef,ngridsx,ngridsy)
end

# TODO: review partitioning of the field
function get_model_loc(parallel_spec::MPISpec)::Tuple{Integer, Integer}
    @unpack ngridsx, ngridsy = parallel_spec
    rank == root_rank && throw("Cannot retrieve model on root rank, you're doing something wrong")
    igrid = (rank%ngridsx)+1
    jgrid = (rank%ngridsy)+1
    return (igrid, jgrid)
end

function update_preconditioner!(model::AbstractModel, ::MPISpec)
    @unpack ngridsx, ngridsy, overlap = model.parallel_spec

    if rank == root_rank
        @info "Full domain for rank $(rank)/$(comm_size)"
    
    else
        igrid = (rank%ngridsx)+1
        jgrid = (rank%ngridsy)+1
        @info "Initialising sub domain model in rank $(rank)/$(comm_size) - $(igrid)x$(jgrid)"

        # TODO: revise grid assignment structure for MPI and BC exchange
        model.parallel_spec.mpiModelArray[igrid,jgrid] = schwarzModel(model;
                                                                      igrid=igrid,
                                                                      jgrid=jgrid,
                                                                      ngridsx=ngridsx,
                                                                      ngridsy=ngridsy, 
                                                                      overlap=overlap)
       
    end
    MPI.Barrier(comm)

    return model
end

function precondition!(model::AbstractModel, ::MPISpec)
    @unpack ngridsx, ngridsy, overlap, niterations, mpiModelArray, damping = model.parallel_spec
    @unpack solver_params = model

    # @info "Preconditioning in rank $(rank)/$(comm_size)"
    if rank == root_rank
        x = WAVI.get_start_guess(model)  
        op = WAVI.get_op(model)
        b = WAVI.get_rhs(model)
        resid = WAVI.get_resid(x,op,b)
        WAVI.set_residual!(model,resid)
        rel_resid = norm(resid)/norm(b)
        converged = rel_resid < solver_params.tol_picard

        model_u_vbuf = UBuffer(model.fields.gu.u, 1)
        model_v_vbuf = UBuffer(model.fields.gv.v, 1)

        MPI.bcast(converged, root_rank, comm)
        model.fields.gu.u[] .= MPI.Scatter!(model_u_vbuf, zeros(size(model.fields.gu.u)), root_rank, comm)
        model.fields.gv.v .= MPI.Scatter!(model_v_vbuf, zeros(size(model.fields.gv.v)), root_rank, comm)
    else
        converged = false
        model_u_vbuf = nothing
        model_v_vbuf = nothing
    end

    if ! converged && rank != root_rank
        igrid, jgrid = get_model_loc(model.parallel_spec)
        model_g = mpiModelArray[igrid, jgrid]

        for iteration = 1:niterations
            schwarzRestrictVelocities!(
                model_g::AbstractModel,
                model::AbstractModel;
                igrid=igrid,
                jgrid=jgrid,
                ngridsx=ngridsx,
                ngridsy=ngridsy,
                overlap=overlap)

            WAVI.update_state!(model_g)

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

    MPI.Barrier(comm)

    if rank == root_rank
        MPI.Gather!(MPI.IN_PLACE, model_u_vbuf, root_rank, comm)
        MPI.Gather!(MPI.IN_PLACE, model_v_vbuf, root_rank, comm)
    end
    return converged, rel_resid
end