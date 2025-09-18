export MPISpec

using JLD2
using MPI
using Parameters
using Plots

using WAVI.Parameters

import WAVI: AbstractField, AbstractGrid, AbstractMeltRate, AbstractModel
import WAVI.Deferred: Collector, register_item!, field_extractor
import WAVI.Fields: GridField, InitialConditions, HGrid, UGrid, VGrid, CGrid, SigmaGrid
import WAVI.Grids: Grid
import WAVI.MeltRates: UniformMeltRate
import WAVI.Models: BasicSpec, Model, get_bed_elevation
import WAVI.Outputs: write_outputs, zip_output, OutputParams
import WAVI.Parameters: TimesteppingParams
import WAVI.Processes: update_state!, update_model_velocities!, update_velocities!, update_velocities_on_h_grid!
import WAVI.Simulations: run_simulation!, timestep!
import WAVI.Time: Clock
import WAVI.Wavelets: UWavelets, VWavelets

# FIXME: important to realise that this specification has become a complex structure to house many things that should be baked into the model structurally
#  not least the global grid and fields 
"""
Struct to represent the MPI parallel specification of a model.

"""
mutable struct MPISpec{N <: Integer, M, G} <: AbstractDecompSpec 
    # MPI Specification information
    px::N
    py::N
    halo::N

    # MPI process information
    global_size::N
    global_comm::M
    global_grid::G
    global_fields::Union{GridField, Nothing}
    field_collector::Collector

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
        @debug "Creating dimensions of $(size) with ($(px), $(py))"
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
            nothing,
            Collector(),
            rank,
            cart_comm,
            [x_coord, y_coord],
            top, right, bottom, left)
    end
end

include("MPI/utils.jl")
include("MPI/exchanges.jl")
include("MPI/outputs.jl")

function Model(grid::G, 
               bed_elevation::Union{Integer, Function, AbstractArray}, 
               spec::S;
               initial_conditions::InitialConditions = InitialConditions(),
               params::Params = Params(),
               solver_params::SolverParams = SolverParams(),
               melt_rate::M = UniformMeltRate()) where {G<:AbstractGrid, S<:MPISpec, M<:AbstractMeltRate}
    @unpack coords, global_size, rank = spec
    @info "[$(rank+1)/$(global_size)] - $(coords) - creating Grid and Model for MPI rank $(rank)"

    # Recalculate grid dimensions and mask parameters, creating a new local Grid
    th, rh, bh, lh = get_halos(spec)
    nx_local, ny_local = get_size(spec)
    x_start, x_end, y_start, y_end = get_bounds(spec)
    
    x0_local = grid.x0 + (x_start-1) * grid.dx
    y0_local = grid.y0 + (y_start-1) * grid.dy
    
    @info "[$(rank+1)/$(global_size)] - proc $(coords[1]),$(coords[2]) - grid $(nx_local)x$(ny_local)"
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

    # dt cannot be copied via the external constructor so we create the structure directly
    local_params = Params(
        params.dt,
        params.g, 
        params.density_ice, 
        params.density_ocean, 
        params.gas_const,
        params.sec_per_year, 
        params.default_thickness, 
        params.default_viscosity,
        params.default_temperature,
        params.default_damage,
        size(params.accumulation_rate) == size(grid)[1:2] ? 
            params.accumulation_rate[x_start:x_end, y_start:y_end] : params.accumulation_rate,
        params.glen_a_activation_energy,
        size(params.glen_a_ref) == size(grid)[1:2] ? 
            params.glen_a_ref[x_start:x_end, y_start:y_end] : params.glen_a_ref,
        params.glen_temperature_ref,
        params.glen_n,
        params.glen_reg_strain_rate,
        size(params.weertman_c) == size(grid)[1:2] ? 
            params.weertman_c[x_start:x_end, y_start:y_end] : params.weertman_c,
        params.weertman_m,
        params.weertman_reg_speed,
        params.sea_level_wrt_geoid,
        params.minimum_thickness,
        params.evolveShelves,
        params.smallHAF
    )

    u_isfixed = grid.u_isfixed[x_start:x_end+1, y_start:y_end]
    v_isfixed = grid.v_isfixed[x_start:x_end, y_start:y_end+1]
    
    # Set halos as fixed velocities
    (spec.left==-1) || (u_isfixed[1,:] .= true; v_isfixed[1,:] .= true)
    (spec.right==-1) || (u_isfixed[end,:] .= true; v_isfixed[end,:] .= true)
    (spec.top==-1) || (u_isfixed[:,1] .= true; v_isfixed[:,1] .= true)
    (spec.bottom==-1) || (u_isfixed[:,end] .= true; v_isfixed[:,end] .= true)

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
    # FIXME: WOAH, bed_array is the wrong size!!! [x_start:x_end, y_start:y_end]
    fields = GridField(local_grid, bed_array[x_start:x_end, y_start:y_end]; initial_conditions=conditions, params=local_params, solver_params)
    model = Model{Float64, Int64, S, GridField, G, M}(local_grid, fields, local_params, solver_params, spec, melt_rate)

    global_bed = typeof(bed_elevation) <: AbstractArray ? bed_elevation : get_bed_elevation(bed_elevation, grid)
    # We provide the full bed as in GridField as it is required for HGrid - this gives us a clean full domain on root
    model.spec.global_fields = GridField(grid, global_bed; initial_conditions, params, solver_params)

    return model
end

##
# Make Model interface work transparently with MPI distributed elements
#

# TODO: remove getproperty, if it's not a global registered field it should be ignored
function Base.getproperty(model::Model{T,N,<:MPISpec,F,G,M}, s::Symbol) where {T,N,F,G,M}
    if s == :global_fields
        # TODO: these need to be registered fields, not user-specified
        ## TODO: fields = collate_global_fields(model.fields, model.spec)
        return model.spec.global_fields
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

##
# Implementations that affect the simuation and data collection
#
function timestep!(model::AbstractModel{T,N,S},
                   timestepping_params::TimesteppingParams,
                   output_params::OutputParams,
                   clock::Clock) where {T,N,S<:MPISpec}
    update_state!(model, clock)

    collect!(model.spec.field_collector, model)

    #write solution if at the first timestep (hack for https://github.com/RJArthern/WAVI.jl/issues/46 until synchronicity is fixed)
    # Have made the interface consistent
    # Have also removed the dependence on individual call
    if (output_params.output_start) && (clock.n_iter == 0)
        write_outputs(model, timestepping_params, output_params, clock)
    end
    
    if timestepping_params.step_thickness
        update_thickness!(model, timestepping_params)
    end
    update_clock!(clock, timestepping_params)

    write_outputs(model, timestepping_params, output_params, clock)

    clear!(model.spec.field_collector)
end

function run_simulation!(model::AbstractModel{T,N,S}, 
                         timestepping_params::TimesteppingParams, 
                         output_params::OutputParams,
                         clock::Clock) where {T,N,S<:MPISpec}
    for field in values(output_params.outputs.items)    
        @info "Registering $(field.path) from outputs"
        if field.path[1] == :global_fields
            register_mpi_field!(model.spec.field_collector, field.path)
        end
    end

    # TODO: we potentially register other fields here too, but currently concentrating on outputs (update_thickness might want to exploit this mechanism)

    for i = (clock.n_iter+1):timestepping_params.n_iter_total
        @info "Running iteration $(clock.n_iter)/$(timestepping_params.n_iter_total)"
        timestep!(model, timestepping_params, output_params, clock)
    end

    zip_output(model, output_params)
end

function register_mpi_field!(collector::Collector, path::Vector{Symbol})
    accessor = function(model)
        collect_mpi_field!(model, path)
        result = model.spec
        for field in path
            result = getproperty(result, field)
        end
        return result
    end
    extractor = field_extractor(join(string.(path), "."), accessor, path)
    register_item!(collector, extractor)
end
