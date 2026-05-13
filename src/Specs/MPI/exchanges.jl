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

function halo_exchange!(model::AbstractModel{<:Any, <:Any, <:MPISpec}; fields::Vector{Symbol}=[:h, :u, :v])
    @unpack halo, rank, comm, top, right, bottom, left = model.spec
    @unpack gh, gu, gv = model.fields

    if halo == 0
        rank == 0 && @warn "No halo exchange to take place, returning"
        return
    end

    th, rh, bh, lh = get_halos(model.spec)

    # Tags for friendly neighbourhood messaging
    top_send_tag = 1
    right_send_tag = 2
    bottom_send_tag = 3
    left_send_tag = 4
    
    # Build field list based on dynamically checking which grid holds the requested fields
    exchange_pairs = Tuple{Any, Symbol}[]
    for f in fields
        hasproperty(gh, f) && push!(exchange_pairs, (gh, f))
        hasproperty(gu, f) && push!(exchange_pairs, (gu, f))
        hasproperty(gv, f) && push!(exchange_pairs, (gv, f))
    end

    # Exchange requested fields
    for (field_data, attribute) in exchange_pairs
        local_field = getproperty(field_data, attribute)

        # We can only perform halo exchange on 2D arrays
        length(size(local_field)) != 2 && continue

        field_nx, field_ny = size(local_field)

        # --- Phase 1: X-Direction Exchange (Left/Right) ---
        requests_x = MPI.RequestSet()

        # Adjust for U-staggering (nx is +1)
        # We need to skip the shared interface face to avoid 1-index shift
        off_x = (field_data === gu) ? 1 : 0

        # Send/Recv Left
        if left > -1
            # Send Left Core: Start at lh+1 + offset
            # Core Face 1 is Interface. Neighbor needs Core Face 2.
            send_left = local_field[lh+1+off_x:lh+halo+off_x, :]
            send_left_flat = copy(reshape(send_left, prod(size(send_left))))
            recv_left_flat = zeros(Float64, prod(size(send_left)))
            push!(requests_x, MPI.Isend(send_left_flat, left, left_send_tag, comm))
            push!(requests_x, MPI.Irecv!(recv_left_flat, left, right_send_tag, comm))
        end

        # Send/Recv Right
        if right > -1
            # Send Right Core: End at field_nx - rh - offset
            # Core Face End is Interface. Neighbor needs Core Face End-1.
            send_right = local_field[field_nx-rh-halo+1-off_x:field_nx-rh-off_x, :]
            send_right_flat = copy(reshape(send_right, prod(size(send_right))))
            recv_right_flat = zeros(Float64, prod(size(send_right)))
            push!(requests_x, MPI.Isend(send_right_flat, right, right_send_tag, comm))
            push!(requests_x, MPI.Irecv!(recv_right_flat, right, left_send_tag, comm))
        end

        MPI.Waitall(requests_x)

        @unpack damping = model.spec

        # Apply X Updates
        if left > -1
            recv_left = reshape(recv_left_flat, halo, field_ny)
            local_field[1:halo, :] .= (1.0 - damping) .* recv_left .+ damping .* local_field[1:halo, :]
        end
        if right > -1
            recv_right = reshape(recv_right_flat, halo, field_ny)
            local_field[field_nx-halo+1:field_nx, :] .= (1.0 - damping) .* recv_right .+ damping .* local_field[field_nx-halo+1:field_nx, :]
        end

        # --- Phase 2: Y-Direction Exchange (Top/Bottom) ---
        requests_y = MPI.RequestSet()

        off_y = (field_data === gv) ? 1 : 0

        # Send/Recv Top
        if top > -1
            # Send Top Core: Start at th+1 + offset
            send_top = local_field[:, th+1+off_y:th+halo+off_y]
            send_top_flat = copy(reshape(send_top, prod(size(send_top))))
            recv_top_flat = zeros(Float64, prod(size(send_top)))
            push!(requests_y, MPI.Isend(send_top_flat, top, top_send_tag, comm))
            push!(requests_y, MPI.Irecv!(recv_top_flat, top, bottom_send_tag, comm))
        end

        # Send/Recv Bottom
        if bottom > -1
            # Send Bottom Core: End at field_ny - bh - offset
            send_bottom = local_field[:, field_ny-bh-halo+1-off_y:field_ny-bh-off_y]
            send_bottom_flat = copy(reshape(send_bottom, prod(size(send_bottom))))
            recv_bottom_flat = zeros(Float64, prod(size(send_bottom)))
            push!(requests_y, MPI.Isend(send_bottom_flat, bottom, bottom_send_tag, comm))
            push!(requests_y, MPI.Irecv!(recv_bottom_flat, bottom, top_send_tag, comm))
        end

        MPI.Waitall(requests_y)

        # Apply Y Updates
        if top > -1
            recv_top = reshape(recv_top_flat, field_nx, halo)
            local_field[:, 1:halo] .= (1.0 - damping) .* recv_top .+ damping .* local_field[:, 1:halo]
        end
        if bottom > -1
            recv_bottom = reshape(recv_bottom_flat, field_nx, halo)
            local_field[:, field_ny-halo+1:field_ny] .= (1.0 - damping) .* recv_bottom .+ damping .* local_field[:, field_ny-halo+1:field_ny]
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

    # Determine global placement for this field's core region (field-aware for staggered grids).
    grid_sym = length(path) >= 2 ? path[2] : :gh
    # Start and end indices for the global field in x
    sx = x_start + lh
    ex = x_end - rh
    # Start and end indices for the global field in y
    sy = y_start + th
    ey = y_end - bh
    # Adjust the end indices for staggered grids
    if grid_sym == :gu
        ex += 1
    elseif grid_sym == :gv
        ey += 1
    elseif grid_sym == :gc
        ex -= 1
        ey -= 1
    end

    # Send/Gather the remote copies from the other nodes into the full field.
    # We provide the local core size and positioning in the target global field.
    field_sz = MPI.Gather(((x_sz - lh - rh, y_sz - th - bh), sx, ex, sy, ey), 0, comm)
    
    if rank == 0
        # We calculate the global grid coordinates for all ranks 
        # based on the received sizes of their core domain (ie. no halo)
        count_sizes = map(x -> prod(x[1]), field_sz)
        field_type = eltype(local_field)
        recv_data = Vector{field_type}(undef, sum(count_sizes))
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
