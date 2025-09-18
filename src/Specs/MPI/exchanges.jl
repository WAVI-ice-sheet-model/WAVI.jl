using Parameters
using Plots
using MPI

using WAVI.Parameters

import WAVI: AbstractField, AbstractGrid, AbstractMeltRate, AbstractModel
import WAVI.Fields: GridField, InitialConditions, HGrid, UGrid, VGrid, CGrid, SigmaGrid
import WAVI.Grids: Grid
import WAVI.MeltRates: UniformMeltRate
import WAVI.Models: BasicSpec, Model, get_bed_elevation
import WAVI.Processes: update_state!, update_model_velocities!, update_velocities!
import WAVI.Wavelets: UWavelets, VWavelets

##
# Additional MPI functionality
#

function halo_exchange!(model::AbstractModel{<:Any, <:Any, <:MPISpec})
    @unpack px, py, halo, global_size, global_comm, rank, comm, coords, top, right, bottom, left = model.spec
    @unpack gh, gu, gv = model.fields
    grid = model.grid

    # @debug "[$(rank+1)/$(global_size)] Halo exchange in progress"
    # Useful to dip out, we'll always be more than one
    if halo == 0
        rank == 0 && @warn "No halo exchange to take place, returning"
        return
    end
    halo_offset = halo - 1          # We need to adjust for use with indexing

    # Tags for neighbourhood messaging
    top_send_tag = 1
    right_send_tag = 2
    bottom_send_tag = 3
    left_send_tag = 4
    
    requests = MPI.RequestSet()
    
    # TODO: currently we hardcode the velocity fields to be exchanged, but this should be user selectable
    for field_sym in [(gu, :u, (1, 0)), (gv, :v, (0, 1))]
        field_data, attribute, add_values = field_sym
        
        local_field = getproperty(field_data, attribute)
        # @debug "[$(rank+1)/$(global_size)] $(size(local_field))"
        x_incr, y_incr = add_values

        # FIXME: POU implementation is under GH#108
        # FIXME: exchange / merging over overlapping fields from both sides of the boundary, POU or otherwise
        
        # Send the "vertical" halo's 
        if left > -1
            send_left = local_field[1:1+halo_offset, :] # FIXME: see above
            send_left_flat = reshape(send_left, prod(size(send_left)))
            recv_left_flat = zeros(Float64, prod(size(send_left)))
            push!(requests, MPI.Isend(send_left_flat, left, left_send_tag, comm))
            push!(requests, MPI.Irecv!(recv_left_flat, left, right_send_tag, comm))
        end

        if right > -1
            send_right = local_field[grid.nx+x_incr-halo_offset:grid.nx+x_incr, :] # FIXME: see above
            send_right_flat = reshape(send_right, prod(size(send_right)))
            recv_right_flat = zeros(Float64, prod(size(send_right)))
            push!(requests, MPI.Isend(send_right_flat, right, right_send_tag, comm))
            push!(requests, MPI.Irecv!(recv_right_flat, right, left_send_tag, comm))
        end

        # Send the "horizontal" halo's
        if top > -1
            send_top = local_field[:, 1:1+halo_offset] # FIXME: see above
            send_top_flat = reshape(send_top, prod(size(send_top)))
            recv_top_flat = zeros(Float64, prod(size(send_top)))
            push!(requests, MPI.Isend(send_top_flat, top, top_send_tag, comm))
            push!(requests, MPI.Irecv!(recv_top_flat, top, bottom_send_tag, comm))
        end

        if bottom > -1
            send_bottom = local_field[:, grid.ny+y_incr-halo_offset:grid.ny+y_incr] # FIXME: see above  
            send_bottom_flat = reshape(send_bottom, prod(size(send_bottom)))
            recv_bottom_flat = zeros(Float64, prod(size(send_bottom)))
            push!(requests, MPI.Isend(send_bottom_flat, bottom, bottom_send_tag, comm))
            push!(requests, MPI.Irecv!(recv_bottom_flat, bottom, top_send_tag, comm))
        end

        # @debug "Waiting on requests"
        MPI.Waitall(requests)
        # @debug "Finished waiting on requests"

        # Update halo regions 
        if left > -1
            recv_left = reshape(recv_left_flat, halo, grid.ny+y_incr)
            local_field[1:1+halo_offset, :] = recv_left # FIXME: see above
        end

        if right > -1
            recv_right = reshape(recv_right_flat, halo, grid.ny+y_incr)
            local_field[grid.nx+x_incr-halo_offset:grid.nx+x_incr, :] = recv_right # FIXME: see above
        end

        if top > -1
            recv_top = reshape(recv_top_flat, grid.nx+x_incr, halo)
            local_field[:, 1:1+halo_offset] = recv_top # FIXME: see above
        end

        if bottom > -1
            recv_bottom = reshape(recv_bottom_flat, grid.nx+x_incr, halo)
            local_field[:, grid.ny+y_incr-halo_offset:grid.ny+y_incr] = recv_bottom # FIXME: see above
        end
    end
end

function collect_mpi_field!(model::AbstractModel{T,N,S}, path::Vector{Symbol}) where {T,N,S<:MPISpec}
    @unpack comm, coords, global_fields, global_size, rank = model.spec

    # Get the full field we want to collect into from the spec, and the equivalent local field on this member
    if path[1] != :global_fields
        error("$(path) should be referring to a global field, so the first symbol should be global_fields")
    end

    global_field = global_fields # Not named correctly
    local_field = model.fields
    for path_el in path[2:end] 
        global_field = getproperty(global_field, path_el)
        local_field = getproperty(local_field, path_el)
    end

    # Establish the local grid information, with full grid information available already from global_grid
    th, rh, bh, lh = get_halos(model.spec)
    # We only handle 2D fields!
    if length(size(local_field)) != 2
        error("Trying to exchange a field ",join(string.(path), ".")," that is not 2D, this is not possible")
    end
    x_sz, y_sz = size(local_field)
    x_start, x_end, y_start, y_end = get_bounds(model.spec)

    @debug "[$(rank+1)/$(global_size)", join(string.(path), "."), "$((x_sz, y_sz, x_start, x_end, y_start, y_end))"

    # Send/Gather the remote copies from the other nodes into the full field 
    # Here we provide the size of the field as well as its positioning in the global grid
    field_sz = MPI.Gather(((x_sz - lh - rh, y_sz - th - bh), 
                           x_start+lh, x_end-rh, y_start+th, y_end-bh), 0, comm)
    
    if rank == 0
        # We calculate the global grid coordinates for all ranks 
        # based on the received sizes of their core domain (ie. no halo)
        count_sizes = map(x -> prod(x[1]), field_sz)
        recv_data = Vector{Float64}(undef, sum(count_sizes))
        recv_buffer = MPI.VBuffer(recv_data, count_sizes)
        @debug "[$(rank+1)/$(global_size) ", join(string.(path), "."), "] Gathering field $((1+lh, size(local_field)[1]-rh, 1+th, size(local_field)[2]-bh)) to buffer $(size(recv_data))"
        MPI.Gatherv!(local_field[1+lh:end-rh, 1+th:end-bh], recv_buffer, comm)

        idxer = collect(cumsum(count_sizes))

        for proc_rank in 0:(global_size-1)
            offset = proc_rank == 0 ? 0 : idxer[proc_rank]
            proc_data = recv_data[offset+1:offset + count_sizes[proc_rank+1]]
            sx, ex, sy, ey = field_sz[proc_rank+1][2:end]
            global_field[sx:ex, sy:ey] = reshape(proc_data, field_sz[proc_rank + 1][1])
        end
    else
        @debug "[$(rank+1)/$(global_size)] Sending ", join(string.(path), "."), " data"
        MPI.Gatherv!(local_field[1+lh:end-rh, 1+th:end-bh], nothing, comm)
    end

    MPI.Barrier(comm)
    return global_field
end
