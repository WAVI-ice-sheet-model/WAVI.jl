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
import WAVI.Processes: update_state!, update_model_velocities!, update_velocities!, update_velocities_on_h_grid!, inner_update!, precondition!, update_rheological_operators!
import WAVI.Simulations: run_simulation!, timestep!
import WAVI.Time: Clock
import WAVI.Wavelets: UWavelets, VWavelets

# FIXME: important to realise that this specification has become a complex structure to house many things that should be baked into the model structurally
#  not least the global grid and fields 
"""
Struct to represent the MPI parallel specification of a model.

Fields:
    px, py: Number of processes in x and y directions
    halo: Number of halo cells on each side of the subgrid
    rank: MPI rank of the current process
    comm: MPI communicator for the current process
    coords: Cartesian coordinates of the current process
    top, right, bottom, left: Neighbouring processes in the x and y directions
    global_size: Total number of processes
    global_comm: MPI communicator for the global grid
    global_grid: Global grid to be used for the model
    damping: Damping factor for halo exchange (default=0.0)
    niterations: Number of Schwarz iterations per Picard iteration (default=5)
    field_collector: Field collector for the model
"""
mutable struct MPISpec{N <: Integer, T <: Number, M <: MPI.Comm, G <: AbstractGrid} <: AbstractDecompSpec
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
    pou::Bool
    damping::T
    niterations::N  # Number of Schwarz iterations per Picard iteration

    @doc """
    Constructor for MPISpec.

    px, py: Number of processes in x and y directions
    halo: Number of halo cells on each side of the subgrid
    grid: Grid to be used for the model
    niterations: Number of Schwarz iterations per Picard iteration (default=5)
    pou: Whether to use partition of unity for halo exchange (default=true)
    damping: Damping factor for halo exchange (default=0.0)
    """
    function MPISpec(px::Integer, py::Integer, halo::Integer, grid::AbstractGrid; pou::Bool=true, damping::AbstractFloat=0.0, niterations::Integer=5)
        (px < 1 || py < 1 || halo < 0) && 
            throw(ArgumentError("Invalid parameters specified for MPISpec"))

        MPI.Initialized() || MPI.Init()
        comm = MPI.COMM_WORLD
        rank = MPI.Comm_rank(comm)
        size = MPI.Comm_size(comm)
        @debug "Creating dimensions of $(size) with ($(px), $(py))"

        if size > 1 && halo == 0
            throw(ArgumentError(
                "MPI Configuration Error: halo must be >= 1 when running with multiple MPI ranks (size=$(size))."
            ))
        end

        if size > 1 && pou && halo < 2
            throw(ArgumentError(
                "MPI Configuration Error: Partition of Unity (pou=true) requires a minimum halo size of 2 for blending."
            ))
        end

        # Create Virtual Cartesian Topology based on no. of procs in each direction
        dims = MPI.Dims_create(size, (px, py))
        cart_comm = MPI.Cart_create(comm, dims)
        x_coord, y_coord = MPI.Cart_coords(cart_comm)

        # Return the rank of the process in the direction specified
        # Note: MPI_Cart_shift returns `MPI_PROC_NULL` (represented as a negative integer)
        #       if the neighbour is not present
        top = MPI.Cart_shift(cart_comm, 1, -1)[2]
        right = MPI.Cart_shift(cart_comm, 0, 1)[2]
        bottom = MPI.Cart_shift(cart_comm, 1, 1)[2]
        left = MPI.Cart_shift(cart_comm, 0, -1)[2]

        @info "[$(rank+1)/$(size)] Neighbours $(top),$(right),$(bottom),$(left)"

        # Check if the process grid is valid
        if px * py != size
            valid = [(i, div(size,i)) for i in 1:size if size % i == 0]
            error(
                "MPI Configuration Error: Invalid process grid.\n\n" *
                "size = $(size), px = $(px), py = $(py)\n\n" *
                "px * py must equal size.\n\n" *
                "Valid (px, py) combinations for size=$(size):\n" *
                join(["  ($(a), $(b))" for (a,b) in valid], "\n")
            )
        end

        # Validation checks for grid size vs halo
        # Need to ensure local grid size is sufficient for halo exchange
        # The constraint is generally: `Core > Halo` (for edge nodes) or `Core > 2*Halo` (for interior nodes)
        # where `Core = div(global_dim, procs)`
        if size != 1
            validate_dimension("nx", grid.nx, px, halo)
            validate_dimension("ny", grid.ny, py, halo)
        end

        return new{Int, Float64, MPI.Comm, AbstractGrid}(
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
            top, right, bottom, left,
            pou,
            damping,
            niterations
            )
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
    @info "[$(rank+1)/$(global_size)] - proc $(coords[1]),$(coords[2]) - $(th), $(rh), $(bh), $(lh)"
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
    # RAS/Schwarz: Only fix the outermost edge cell - this provides Dirichlet BC from neighbour
    # Interior halo cells (2:halo) are SOLVED locally, then DISCARDED during exchange
    # The simple overwrite in halo_exchange! replaces them with neighbour's core values
    (spec.left < 0) || (u_isfixed[1,:] .= true; v_isfixed[1,:] .= true)
    (spec.right < 0) || (u_isfixed[end,:] .= true; v_isfixed[end,:] .= true)
    (spec.top < 0) || (u_isfixed[:,1] .= true; v_isfixed[:,1] .= true)
    (spec.bottom < 0) || (u_isfixed[:,end] .= true; v_isfixed[:,end] .= true)

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

    if typeof(bed_elevation) <: AbstractArray
        bed_array = bed_elevation[x_start:x_end, y_start:y_end]
    else
        bed_array = get_bed_elevation(bed_elevation, local_grid)
    end

    # Create local mpi_rank field filled with this rank's number
    local_mpi_rank = fill(Float64(rank), nx_local, ny_local)

    fields = GridField(local_grid, bed_array; initial_conditions=conditions, params=local_params, solver_params, mpi_rank=local_mpi_rank)
    model = Model{Float64, Int64, S, GridField, G, M}(local_grid, fields, local_params, solver_params, spec, melt_rate)

    global_bed = typeof(bed_elevation) <: AbstractArray ? bed_elevation : get_bed_elevation(bed_elevation, grid)
    # Create global mpi_rank field (will be populated during collection)
    global_mpi_rank = zeros(Float64, grid.nx, grid.ny)

    # We provide the full bed as in GridField as it is required for HGrid - this gives us a clean full domain on root
    model.spec.global_fields = GridField(grid, global_bed; initial_conditions, params, solver_params, mpi_rank=global_mpi_rank)

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
    update_velocities!(model)
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

    #write solution if at the first timestep (hack for https://github.com/RJArthern/WAVI.jl/issues/46 until synchronicity is fixed)
    # Have made the interface consistent
    # Have also removed the dependence on individual call
    if (output_params.output_start) && (clock.n_iter == 0)
        collect!(model.spec.field_collector, model)
        write_outputs(model, timestepping_params, output_params, clock)
        clear!(model.spec.field_collector)
    end
    
    if timestepping_params.step_thickness
        update_thickness!(model, timestepping_params)
        # Sync thickness halos; do not RAS-overwrite u/v when using PoU (see mpi_sync_halos_after_thickness!)
        mpi_sync_halos_after_thickness!(model)
    end
    update_clock!(clock, timestepping_params)

    # Collect AFTER thickness update to match BasicSpec output timing
    collect!(model.spec.field_collector, model)
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

    mpi_sync_halos_initial!(model)

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

"""
    inner_update!(model::Model{<:Any, <:Any, <:MPISpec})

Overload for the inner update of the velocity solve.
We need to sync halos before performing the update so that the viscosity and other
calculations have correct boundary information from neighbouring procs.

Note: `update_rheological_operators!` is called a second time here (after the halo sync)
because the base `inner_update!` builds the operators before we have correct cross-rank
rheology values (β, βeff, ηav). The second call rebuilds them with consistent halo data.
"""
function inner_update!(model::Model{<:Any, <:Any, <:MPISpec})
    # RAS ghost sync for velocities; AS-PoU keeps overlap values from the last prolongation.
    if !model.spec.pou
        halo_exchange!(model; fields=[:u, :v])
    end
    # Call the standard inner update function for the velocity solve
    invoke(inner_update!, Tuple{AbstractModel}, model)
    # Sync rheology fields used near rank interfaces, then rebuild operators with correct halo data.
    halo_exchange!(model; fields=[:β, :βeff, :ηav])
    update_rheological_operators!(model)
    return model
end

function core_inner_masks(model::Model{<:Any, <:Any, <:MPISpec})
    @unpack gu, gv = model.fields
    th, rh, bh, lh = get_halos(model.spec)

    u_core_mask = falses(size(gu.mask_inner))
    u_core_mask[(1+lh):(size(u_core_mask, 1)-rh), (1+th):(size(u_core_mask, 2)-bh)] .= true
    u_core_inner = u_core_mask[gu.mask_inner]

    v_core_mask = falses(size(gv.mask_inner))
    v_core_mask[(1+lh):(size(v_core_mask, 1)-rh), (1+th):(size(v_core_mask, 2)-bh)] .= true
    v_core_inner = v_core_mask[gv.mask_inner]

    return u_core_inner, v_core_inner
end

"""
    precondition!(model::Model{<:Any, <:Any, <:MPISpec})

Solves the linear system using an iterative overlapping Schwarz method across the distributed domain.

**Iterate** up to `niterations` times:
*   **Local Solve**: Applies the standard `precondition!` locally on each rank.
*   **Interface Update**: Exchanges updated velocities with neighbouring ranks. 
    Depending on the model configuration, this uses either Additive Schwarz with Partition-of-Unity 
    (`pou = true`) or standard Restricted Additive Schwarz (`pou = false`).
*   **Check Convergence**: Computes the global relative residual at the end of each iteration.
    To avoid double-counting, only the non-overlapping "core" unknowns are used in the global MPI reduction.
    *   Each rank computes the squared norm of its local residual contribution using **core unknowns only**
        (excluding overlap/halo regions) to prevent double-counting.
    *   Values are aggregated across all processes using `MPI.Allreduce` with `MPI.SUM`.
    *   Solver exits early if the global relative residual meets the Picard tolerance.
"""
function precondition!(model::Model{<:Any, <:Any, <:MPISpec})
    @unpack niterations, pou = model.spec
    @unpack solver_params = model

    converged = false
    global_rel_resid = Inf

    for iteration = 1:niterations
        if (iteration > 1) && (model.spec.rank == 0)
            @debug "Schwarz iteration $iteration"
        end

        if pou
            # Snapshot current velocities before the local solve.
            # These are used later to apply Additive Schwarz damping and 
            # Partition-of-Unity (PoU) weighting during the interface update.
            u0 = copy(model.fields.gu.u)
            v0 = copy(model.fields.gv.v)
            invoke(precondition!, Tuple{AbstractModel}, model)
            mpi_pou_weighted_prolong_velocities!(model, u0, v0)
        else
            invoke(precondition!, Tuple{AbstractModel}, model)
            halo_exchange!(model; fields=[:u, :v])
        end

        # Check convergence after each iteration
        x = get_start_guess(model)
        op = get_op(model)
        b = get_rhs(model)
        resid = get_resid(x, op, b)

        # Global Residual Check (core-only):
        # exclude overlap halos so each physical unknown is counted once globally.
        @unpack gu, gv = model.fields
        u_core_inner, v_core_inner = core_inner_masks(model)

        # Calculate squared norms locally on core unknowns only
        local_resid_sq = sum(abs2, @view resid[1:gu.ni][u_core_inner]) +
                         sum(abs2, @view resid[(gu.ni+1):(gu.ni+gv.ni)][v_core_inner])
        local_rhs_sq = sum(abs2, @view b[1:gu.ni][u_core_inner]) +
                       sum(abs2, @view b[(gu.ni+1):(gu.ni+gv.ni)][v_core_inner])

        # Sum squared norms across all ranks
        global_resid_sq = MPI.Allreduce(local_resid_sq, MPI.SUM, model.spec.comm)
        global_rhs_sq = MPI.Allreduce(local_rhs_sq, MPI.SUM, model.spec.comm)

        # Compute global relative residual (guard against degenerate zero-RHS case)
        global_rel_resid = iszero(global_rhs_sq) ? sqrt(global_resid_sq) : sqrt(global_resid_sq) / sqrt(global_rhs_sq)

        converged = global_rel_resid < solver_params.tol_picard
        if converged
            if model.spec.rank == 0
                @debug "Schwarz converged early at iteration $iteration"
            end
            break
        end
    end

    @info "Picard Check: Global Relative Residual = $global_rel_resid (Tol = $(solver_params.tol_picard))"

    return converged, global_rel_resid
end
