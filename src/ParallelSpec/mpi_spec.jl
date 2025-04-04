export MPISpec, update_preconditioner!, precondition!

"""
Naive implementation of a workflow for coordination of domain decomposition solving 
across MPI workers

See GH#91

"""

using MPI
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm) + 1
comm_size = MPI.Comm_size(comm)
root_rank = 1

@with_kw struct MPISpec{T <: Real, N <: Integer} <: AbstractParallelSpec
    ngridsx::N = 1 
    ngridsy::N = 1
    overlap::N = 0
    niterations::N = 0
    damping::T = 0.0
    mpiModelArray::Array{AbstractModel,2} = Array{AbstractModel,2}(undef,ngridsx,ngridsy)
end

"""
update_velocities! for MPISpec

Solve momentum equation to update the velocities on rank, coordinate BC transfer from root

"""
function update_velocities!(model::AbstractModel{T,N,M,PS}) where {T<:Real, N<:Integer, M<:AbstractMeltRate, PS<:MPISpec}
    @unpack params,solver_params=model
    @unpack gu,gv,wu,wv = model.fields

    # 1. Prepare all "mini-WAVI" instances across the domain, initialise our global domain representation
    if (rank == root_rank)
        @info "update_veloctities! for MPI, adapting conditioner"
        update_preconditioner!(model)
    end

    MPI.barrier()

    # 2. Solve independently across the subdomains and send back to the root domain
    # TODO: refactor for communication of 
    # TODO: incorporate transfer of halos between iterations
    if (rank != root_rank)
        converged::Bool = false
        i_picard::Int64 = 0
        rel_resid = Inf

        while !converged && (i_picard < solver_params.maxiter_picard)
            i_picard = i_picard + 1
            inner_update!(model)
            converged, rel_resid = precondition!(model)
        end

        println("Solved momentum equation on rank $(rank)/$(comm_size) with residual ", 
            round(rel_resid, sigdigits=3)," at iteration ", i_picard)
    end

    MPI.barrier()

    if (rank == root_rank)

    end
    
    return model
end


function update_preconditioner!(model::AbstractModel, ::MPISpec)
    @unpack ngridsx, ngridsy, overlap = model.parallel_spec

    if rank == root_rank
        @info "Full domain for rank $(rank)/$(comm_size)"
    
    else
        igrid = (rank%ngridsx)+1
        jgrid = (rank%ngridsy)+1
        @info "Initialising sub domain model in rank $(rank)/$(comm_size) - $(igrid)x$(jgrid)"

        model.parallel_spec.mpiModelArray[igrid,jgrid] = schwarzModel(model;
                                                                      igrid=igrid,
                                                                      jgrid=jgrid,
                                                                      ngridsx=ngridsx,
                                                                      ngridsy=ngridsy, 
                                                                      overlap=overlap)
       
    end

    return model
end

function precondition!(model::AbstractModel, ::MPISpec)

end