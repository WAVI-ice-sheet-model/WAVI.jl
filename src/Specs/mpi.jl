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

include("MPI/exchanges.jl")
include("MPI/utils.jl")

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

#Base.propertynames(model::Model{T,N,<:MPISpec,F,G,M}, private::Bool) where {T,N,F,G,M} = (fieldnames(typeof(model)..., :global_fields))
Base.propertynames(model::Model{T,N,S,F,G,M}, private::Bool) where {T,N,S,F,G,M} = (fieldnames(typeof(model)..., :global_fields))

# FIXME: this is a sign of a frustration in WAVIs structural layout - too many deep nested structures accessed through high level passing
#  which inhibits multiple dispatch
function Base.getproperty(model::Model{T,N,S,F,G,M}, s::Symbol) where {T,N,S,F,G,M}
#    @info "getproperty from generic model spec"
    if s == :global_fields
#        @info "Calling global fields on local model"
        return getfield(model, :fields)
    elseif s == :global_grid
#        @info "Calling global grid on local model"
        return getfield(model, :grid)
    end
    getfield(model, s)
end

function Base.getproperty(model::Model{T,N,<:MPISpec,F,G,M}, s::Symbol) where {T,N,F,G,M}
#    @info "getproperty from MPI model spec $(s)"
    if s == :global_fields
        # TODO: these need to be registered fields, not user-specified
        fields = collate_global_fields(model.fields, model.spec)
        return fields
    elseif s == :global_grid
#        @info "Using global grid from original MPI spec"
        return model.spec.global_grid
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



