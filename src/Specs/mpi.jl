export MPISpec

using JLD2
using Parameters
using Plots
using MPI

using WAVI.Parameters

import WAVI: AbstractField, AbstractGrid, AbstractMeltRate, AbstractModel
import WAVI.Fields: GridField, InitialConditions, HGrid, UGrid, VGrid, CGrid, SigmaGrid
import WAVI.Grids: Grid
import WAVI.MeltRates: UniformMeltRate
import WAVI.Models: BasicSpec, Model, get_bed_elevation
import WAVI.Outputs: write_outputs, zip_output, OutputParams
import WAVI.Parameters: TimesteppingParams
import WAVI.Processes: update_state!, update_model_velocities!, update_velocities!
import WAVI.Time: Clock
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

# TODO: remove getproperty, if it's not a global registered field it should be ignored
function Base.getproperty(model::Model{T,N,<:MPISpec,F,G,M}, s::Symbol) where {T,N,F,G,M}
    if s == :global_fields
        # TODO: these need to be registered fields, not user-specified
        fields = collate_global_fields(model.fields, model.spec)
        return fields
    elseif s == :global_grid
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

# TODO: these redefinitions are loathsome, but can't get a clearer way of limiting output via a simple method declaration with @root
function write_output(model::M, output_params::OutputParams, clock::Clock) where {M<:AbstractModel{<:Any, <:Any, <:MPISpec}}
    output_dict = collect!(output_params, model)
    name = lpad(clock.n_iter, 10,"0")

    @root begin
        if isnothing(output_dict)
            @warn "No outputs processed for $(name)"
        end

        #put the grid co-ordinates and time into output.
        #Round time in output to some decimal places to make it prettier (machine precision can make this look nasty!)
        if ~haskey(output_dict, :t); output_dict["t"] = round(clock.time, digits = 3); end
        if ~haskey(output_dict, :x); output_dict["x"] = model.global_grid.xxh; end
        if ~haskey(output_dict, :y); output_dict["y"] = model.global_grid.yyh; end

        fname = string(output_params.output_path, output_params.prefix, name)
        if output_params.output_format == "jld2"
            fname = string(fname, ".jld2")
            save(fname, output_dict)
        elseif output_params.output_format == "mat"
            fname = string(fname, ".mat")
            matwrite(fname, output_dict)
        end
        
        @info "Output at timestep number $(clock.n_iter) - $(fname)"
    end
    clear!(output_params)
end

function write_outputs(model::M,
                       timestepping_params::TimesteppingParams, 
                       output_params::OutputParams, 
                       clock::Clock) where {M<:AbstractModel{<:Any, <:Any, <:MPISpec}}
    @root begin
        #check if we have hit a permanent checkpoint
        if mod(clock.n_iter, timestepping_params.n_iter_chkpt) == 0
            #output a permanent checkpoint
            n_iter_string =  lpad(clock.n_iter, 10, "0"); #filename as a string with 10 digits
            fname = joinpath(output_params.output_path, string("Chkpt_",n_iter_string, ".jld2"))
            @save fname model=model timestepping_params=timestepping_params clock=clock
            @info "MPI permanent checkpoint at timestep number $(clock.n_iter)"
        end
    end

    #check if we have hit an output timestep
    if mod(clock.n_iter, output_params.n_iter_out) == 0
        write_output(model, output_params, clock)
    end

        #check the dump velocity flag at the final timestep
#         if (clock.n_iter == timestepping_params.n_iter_total) && output_params.dump_vel
#             write_vel(output_params, model)
#         end
end

function zip_output(model::M, output_params::OutputParams) where {M<:AbstractModel{<:Any, <:Any, <:MPISpec}}
    @info "[$(model.spec.rank + 1)/$(model.spec.global_size)] Called zip output"
    @root begin
        if output_params.zip_format == "nc"
            nc_name_full = string(output_params.output_path, output_params.prefix, ".nc")
            @info "Creating NetCDF output on MPI root $(nc_name_full)"
            make_ncfile(output_params.output_format, output_params.output_path, nc_name_full, output_params.prefix)
        end
    end
    return nothing
end

