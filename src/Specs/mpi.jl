export MPISpec

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

include("MPI/utils.jl")

struct MPISpec{N <: Integer, M, G} <: AbstractDecompSpec 
    # MPI Specification information
    px::N
    py::N
    halo::N

    # MPI process information
    global_size::N
    global_comm::M
    global_grid::G

    rank::N
    comm::M
    coords::Array{N, 1}

    # Neighbourhood information
    top::Union{N, Nothing}
    right::Union{N, Nothing}
    bottom::Union{N, Nothing}
    left::Union{N, Nothing}

    function MPISpec(px::Int, py::Int, halo::Int, grid::AbstractGrid)
        (px < 1 || py < 1 || halo < 0) && 
            throw(ArgumentError("Invalid parameters specified for MPISpec"))

        MPI.Init()
        comm = MPI.COMM_WORLD
        rank = MPI.Comm_rank(comm)
        size = MPI.Comm_size(comm)
        
        dims = MPI.Dims_create(size, (px, py))
        cart_comm = MPI.Cart_create(comm, dims)
        x_coord, y_coord = MPI.Cart_coords(cart_comm)

        top = MPI.Cart_shift(cart_comm, 1, -1)[2]
        right = MPI.Cart_shift(cart_comm, 0, 1)[2]
        bottom = MPI.Cart_shift(cart_comm, 1, 1)[2]
        left = MPI.Cart_shift(cart_comm, 0, -1)[2]

        @info "[$(rank+1)/$(size)] Neighbours $(top),$(right),$(bottom),$(left)"
        
        return new{Int, MPI.Comm, AbstractGrid}(
            px, 
            py, 
            halo,
            size,
            comm,
            grid,
            rank,
            cart_comm,
            [x_coord, y_coord],
            top, right, bottom, left)
    end
end

function Model(grid::G, 
               bed_elevation::Union{Integer, Function, AbstractArray}, 
               spec::S;
               initial_conditions::InitialConditions = InitialConditions(),
               params::Params = Params(),
               solver_params::SolverParams = SolverParams(),
               melt_rate::M = UniformMeltRate()) where {G<:AbstractGrid, S<:MPISpec, M<:AbstractMeltRate}
    @unpack px, py, halo, global_size, global_comm, rank, comm, coords, top, right, bottom, left = spec

    @info "[$(rank+1)/$(global_size)] - $(coords) - creating Grid and Model for MPI rank $(rank)"
    x_coord, y_coord = coords

    # Recalculate grid dimensions and mask parameters, creating a new local Grid
    th, rh, bh, lh = top > -1 ? halo : 0, right > -1 ? halo : 0, bottom > -1 ? halo : 0, left > -1 ? halo : 0
    nx_local = div(grid.nx, px) + lh + rh
    ny_local = div(grid.ny, py) + th + bh
    
    x_start = max(x_coord * div(grid.nx, px) + 1 - lh, 1)
    y_start = max(y_coord * div(grid.ny, py) + 1 - th, 1)

    x_end = min(x_coord * div(grid.nx, px) + div(grid.nx, px) + rh, grid.nx)
    y_end = min(y_coord * div(grid.ny, py) + div(grid.ny, py) + bh, grid.ny)
    
    x0_local = grid.x0 + (x_start-1) * grid.dx
    y0_local = grid.y0 + (y_start-1) * grid.dy
    
    @info "[$(rank+1)/$(global_size)] - proc $(x_coord),$(y_coord) - grid $(nx_local)x$(ny_local)"
    @info "[$(rank+1)/$(global_size)] - X [$(x_start):$(x_end)] - Y [$(y_start):$(y_end)] - Centroid $(x0_local),$(y0_local) "

    u_grid_size, v_grid_size = (grid.nx+1, grid.ny), (grid.nx, grid.ny+1)
    
    conditions = InitialConditions(
        initial_thickness = size(initial_conditions.initial_thickness) == size(grid)[1:2] ? 
            initial_conditions.initial_thickness[x_start:x_end, y_start:y_end] : initial_conditions.initial_thickness,
        initial_grounded_fraction = size(initial_conditions.initial_grounded_fraction) == size(grid)[1:2] ? 
            initial_conditions.initial_grounded_fraction[x_start:x_end, y_start:y_end] : initial_conditions.initial_grounded_fraction,
        initial_u_veloc = size(initial_conditions.initial_u_veloc) == u_grid_size ? 
            initial_conditions.initial_u_veloc[x_start:x_end+1, y_start:y_end] : initial_conditions.initial_u_veloc,
        initial_v_veloc = size(initial_conditions.initial_v_veloc) == v_grid_size ? 
            initial_conditions.initial_v_veloc[x_start:x_end, y_start:y_end+1] : initial_conditions.initial_v_veloc,
        initial_viscosity = size(initial_conditions.initial_viscosity) == size(grid) ? 
            initial_conditions.initial_viscosity[x_start:x_end, y_start:y_end, :] : initial_conditions.initial_viscosity,
        initial_temperature = size(initial_conditions.initial_temperature) == size(grid) ? 
            initial_conditions.initial_temperature[x_start:x_end, y_start:y_end, :] : initial_conditions.initial_temperature,
        initial_damage = size(initial_conditions.initial_damage) == size(grid) ? 
            initial_conditions.initial_damage[x_start:x_end, y_start:y_end, :] : initial_conditions.initial_damage
    )

    u_isfixed = grid.u_isfixed[x_start:x_end+1, y_start:y_end]
    v_isfixed = grid.v_isfixed[x_start:x_end, y_start:y_end+1]
    
    # Set halos as fixed velocities
    (left==-1) || (u_isfixed[1:1+lh,:] .= true; v_isfixed[1:1+lh,:] .= true)
    (right==-1) || (u_isfixed[end-rh:end,:] .= true; v_isfixed[end-rh:end,:] .= true)
    (top==-1) || (u_isfixed[:,1:1+th] .= true; v_isfixed[:,1:1+th] .= true)
    (bottom==-1) || (u_isfixed[:,end:end-bh] .= true; v_isfixed[:,end:end-bh] .= true)

    local_grid = Grid(
        nx = nx_local,
        ny = ny_local,
        nσ = grid.nσ,
        dx = grid.dx,
        dy = grid.dy,
        x0 = x0_local,
        y0 = y0_local,
        h_mask = grid.h_mask[x_start:x_end, y_start:y_end],
        h_isfixed = grid.h_isfixed[x_start:x_end, y_start:y_end],
        u_iszero = grid.u_iszero[x_start:x_end+1, y_start:y_end],
        v_iszero = grid.v_iszero[x_start:x_end, y_start:y_end+1],
        u_isfixed = u_isfixed,
        v_isfixed = v_isfixed,
        quadrature_weights = grid.quadrature_weights,
        σ = grid.σ,
        basin_ID = grid.basin_ID[x_start:x_end, y_start:y_end])

    bed_array = get_bed_elevation(bed_elevation, local_grid)
    fields = GridField(local_grid, bed_array; initial_conditions=conditions, params, solver_params)
    model = Model{Float64, Int64, S, GridField, G, M}(local_grid, fields, params, solver_params, spec, melt_rate)
    return model
end

##
# Make Model interface work transparently with MPI distributed elements
#

Base.propertynames(model::Model{T,N,<:MPISpec,F,G,M}, private::Bool) where {T,N,F,G,M} = (fieldnames(typeof(model)..., :global_fields))

function Base.getproperty(model::Model{T,N,<:MPISpec,F,G,M}, s::Symbol) where {T,N,F,G,M}
    if s == :global_fields
        fields = collate_global_fields(model.fields, model.spec)
        return fields
    end
    return getfield(model, s)
end

##
# Override Model oriented methods to intercept calls that need extra processing
#

function update_model_velocities!(model::Model{<:Any, <:Any, <:MPISpec})
    @unpack px, py, halo, global_size, global_comm, rank, comm, coords = model.spec
    update_velocities!(model)

    @debug "[$(rank+1)/$(global_size)] - hitting velocity solve barrier, exchanging halos"
    MPI.Barrier(comm)

    halo_exchange!(model)
    return model
end


##
# Overrides for other functionalities that need restricting to root node
#
# TODO: override @debug, @info, @warn and @error for MPI based logging, with the rank out of size and / or grid location



##
# Additional MPI functionality
#

function halo_exchange!(model::Model{<:Any, <:Any, <:MPISpec})
    @unpack px, py, halo, global_size, global_comm, rank, comm, coords, top, right, bottom, left = model.spec
    @unpack gh, gu, gv = model.fields
    grid = model.grid

    @debug "[$(rank+1)/$(global_size)] Halo exchange in progress"
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
    
    for field_sym in [(gu, :u, (1, 0)), (gv, :v, (0, 1))]
        field_data, attribute, add_values = field_sym
        
        local_field = getproperty(field_data, attribute)
        # @debug "[$(rank+1)/$(global_size)] $(size(local_field))"
        x_incr, y_incr = add_values

        # TODO: got a headache with the overloading of naming conventions - recheck all this diagrammatically and conform to WAVI
        # TODO: need to make sure halos are correctly configured - corners go to both neighbours is it?

        # Send the "vertical" halo's 
        if left > -1
            @debug "[$(rank+1)/$(global_size)] Sending left to $(left)"
            send_left = local_field[1:1+halo_offset, :]
            send_left_flat = reshape(send_left, prod(size(send_left)))
            recv_left_flat = zeros(Float64, prod(size(send_left)))
            push!(requests, MPI.Isend(send_left_flat, left, left_send_tag, comm))
            push!(requests, MPI.Irecv!(recv_left_flat, left, right_send_tag, comm))
        end

        if right > -1
            @debug "[$(rank+1)/$(global_size)] Sending right to $(right)"
            send_right = local_field[grid.nx+x_incr-halo_offset:grid.nx+x_incr, :]
            send_right_flat = reshape(send_right, prod(size(send_right)))
            recv_right_flat = zeros(Float64, prod(size(send_right)))
            push!(requests, MPI.Isend(send_right_flat, right, right_send_tag, comm))
            push!(requests, MPI.Irecv!(recv_right_flat, right, left_send_tag, comm))
        end

        # Send the "horizontal" halo's
        if top > -1
            @debug "[$(rank+1)/$(global_size)] Sending top to $(top)"
            send_top = local_field[:, 1:1+halo_offset]
            send_top_flat = reshape(send_top, prod(size(send_top)))
            recv_top_flat = zeros(Float64, prod(size(send_top)))
            push!(requests, MPI.Isend(send_top_flat, top, top_send_tag, comm))
            push!(requests, MPI.Irecv!(recv_top_flat, top, bottom_send_tag, comm))
        end

        if bottom > -1
            @debug "[$(rank+1)/$(global_size)] Sending bottom to $(bottom)"
            send_bottom = local_field[:, grid.ny+y_incr-halo_offset:grid.ny+y_incr]        
            send_bottom_flat = reshape(send_bottom, prod(size(send_bottom)))
            recv_bottom_flat = zeros(Float64, prod(size(send_bottom)))
            push!(requests, MPI.Isend(send_bottom_flat, bottom, bottom_send_tag, comm))
            push!(requests, MPI.Irecv!(recv_bottom_flat, bottom, top_send_tag, comm))
        end

        @debug "Waiting on requests"
        MPI.Waitall(requests)
        @debug "Finished waiting on requests"

        # Update halo regions 
        if left > -1
            recv_left = reshape(recv_left_flat, halo, grid.ny+y_incr)
            local_field[1:1+halo_offset, :] = recv_left
        end

        if right > -1
            recv_right = reshape(recv_right_flat, halo, grid.ny+y_incr)
            local_field[grid.nx+x_incr-halo_offset:grid.nx+x_incr, :] = recv_right
        end

        if top > -1
            recv_top = reshape(recv_top_flat, grid.nx+x_incr, halo)
            local_field[:, 1:1+halo_offset] = recv_top
        end

        if bottom > -1
            recv_bottom = reshape(recv_bottom_flat, grid.nx+x_incr, halo)
            local_field[:, grid.ny+y_incr-halo_offset:grid.ny+y_incr] = recv_bottom
        end
    end
end

function collate_global_fields(local_fields::AbstractField, spec::MPISpec)
    @unpack px, py, halo, global_size, global_comm, rank, comm, coords, top, right, bottom, left, global_grid = spec

    grid = global_grid
    x_coord, y_coord = coords

    @debug "[$(rank+1)/$(global_size)] Generating local version of the global grid for collation at $(grid.nx) x $(grid.ny)"

    gh = HGrid(nxh=grid.nx, nyh=grid.ny, b=zeros(grid.nx, grid.ny))
    gu = UGrid(nxu=grid.nx+1, nyu=grid.ny, dx=grid.dx, dy=grid.dy, levels=local_fields.gu.levels)
    gv = VGrid(nxv=grid.nx, nyv=grid.ny+1, dx=grid.dx, dy=grid.dy, levels=local_fields.gv.levels)
    gc = CGrid(nxc=grid.nx-1, nyc=grid.ny-1)

    g3d = SigmaGrid(
        nxs = grid.nx,
        nys = grid.ny,
        nσs = grid.nσ,
        σ = grid.σ,
        η = zeros(grid.nx, grid.ny, grid.nσ),
        θ = zeros(grid.nx, grid.ny, grid.nσ),
        Φ = zeros(grid.nx, grid.ny, grid.nσ),
        glen_b = zeros(grid.nx, grid.ny, grid.nσ),
        quadrature_weights = grid.quadrature_weights
    )    
    # Wavelet-grids
    wu = UWavelets(nxuw=grid.nx+1, nyuw=grid.ny, levels=local_fields.gu.levels)
    wv = VWavelets(nxvw=grid.nx, nyvw=grid.ny+1, levels=local_fields.gv.levels)

    ## TODO: duplicated from above, that's not a good smell
    th, rh, bh, lh = top > -1 ? halo : 0, right > -1 ? halo : 0, bottom > -1 ? halo : 0, left > -1 ? halo : 0
    x_start = max(x_coord * div(grid.nx, px) + 1 - lh, 1)
    y_start = max(y_coord * div(grid.ny, py) + 1 - th, 1)

    x_end = min(x_coord * div(grid.nx, px) + div(grid.nx, px) + rh, grid.nx)
    y_end = min(y_coord * div(grid.ny, py) + div(grid.ny, py) + bh, grid.ny)
    # end of clone

    #Wavelet-grid, u-component.
    wu=UWavelets(nxuw=grid.nx+1,nyuw=grid.ny,levels=local_fields.wu.levels)

    #Wavelet-grid, v-component.
    wv=VWavelets(nxvw=grid.nx,nyvw=grid.ny+1,levels=local_fields.wv.levels)
    global_fields = GridField(gh,gu,gv,gc,g3d,wu,wv)

    MPI.Barrier(comm)  
    @debug "[$(rank+1)/$(global_size)] - LOCAL U $(size(local_fields.gu.u)) GLOBAL [$(x_start):$(x_end+1), $(y_start):$(y_end)]"
    @debug "[$(rank+1)/$(global_size)] - LOCAL V $(size(local_fields.gv.v)) GLOBAL [$(x_start):$(x_end), $(y_start):$(y_end+1)]"
    @debug "[$(rank+1)/$(global_size)] Commencing data transfers"
  
    x_sz, y_sz = size(local_fields.gu.u)
    
    # Remove halos and ensure that, if not at the edge of the global boundary, we're not including the extra row
    u_grids = MPI.Gather(((x_sz - lh - rh, y_sz - th - bh), 
                         x_start+lh, x_end-rh+1, y_start+th, y_end-bh), 0, comm)

    if rank == 0
        # We calculate the global grid coordinates for all ranks 
        # based on the received sizes of their core domain (ie. no halo)
        count_sizes = map(x -> prod(x[1]), u_grids)
        recv_data = Vector{Float64}(undef, sum(count_sizes))
        recv_buffer = MPI.VBuffer(recv_data, count_sizes)
        MPI.Gatherv!(local_fields.gu.u[1+lh:end-rh, 1+th:end-bh], recv_buffer, comm)
        
        idxer = collect(cumsum(count_sizes))

        for proc_rank in 0:(global_size-1)
            offset = proc_rank == 0 ? 0 : idxer[proc_rank]
            proc_data = recv_data[offset+1:offset + count_sizes[proc_rank+1]]
            sx, ex, sy, ey = u_grids[proc_rank+1][2:end]
            global_fields.gu.u[sx:ex, sy:ey] = reshape(proc_data, u_grids[proc_rank + 1][1])
        end
    else
        MPI.Gatherv!(local_fields.gu.u[1+lh:end-rh, 1+th:end-bh], nothing, comm)
    end

    MPI.Barrier(comm)

    x_sz, y_sz = size(local_fields.gv.v)
    v_grids = MPI.Gather(((x_sz - lh - rh, y_sz - th - bh), 
                         x_start+lh, x_end-rh, y_start+th, y_end-bh+1), 0, comm)
    
    if rank == 0
        # We calculate the global grid coordinates for all ranks 
        # based on the received sizes of their core domain (ie. no halo)
        count_sizes = map(x -> prod(x[1]), v_grids)
        recv_data = Vector{Float64}(undef, sum(count_sizes))
        recv_buffer = MPI.VBuffer(recv_data, count_sizes)
        MPI.Gatherv!(local_fields.gv.v[1+lh:end-rh, 1+th:end-bh], recv_buffer, comm)
        
        idxer = collect(cumsum(count_sizes))

        for proc_rank in 0:(global_size-1)
            offset = proc_rank == 0 ? 0 : idxer[proc_rank]
            proc_data = recv_data[offset+1:offset + count_sizes[proc_rank+1]]
            sx, ex, sy, ey = v_grids[proc_rank+1][2:end]
            global_fields.gv.v[sx:ex, sy:ey] = reshape(proc_data, v_grids[proc_rank + 1][1])
        end
    else
        MPI.Gatherv!(local_fields.gv.v[1+lh:end-rh, 1+th:end-bh], nothing, comm)
    end

    MPI.Barrier(comm)
    return global_fields
end