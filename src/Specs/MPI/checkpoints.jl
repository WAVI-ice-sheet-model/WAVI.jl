# MPISpec checkpointing: gather fields to root, save one serial BasicSpec file,
# and on pickup scatter that file into the live MPI model.

using WAVI.Outputs: write_checkpoint!, load_checkpoint
import WAVI.Fields: initial_conditions_from_fields
import WAVI.Processes: update_rheological_operators!, update_velocities_on_h_grid!,
    update_model_wavelets!, update_shelf_strain_rate!, update_surf_speed!

# Paths into a GridField (works for both local `fields` and `global_fields`).
# Each path is a tuple of property names, e.g. (:gh, :h) means fields.gh.h.
const CHECKPOINT_INITIAL_CONDITION_FIELD_PATHS = (
    (:gh, :h),
    (:gh, :grounded_fraction),
    (:gu, :u),
    (:gv, :v),
    (:g3d, :η),
    (:g3d, :θ),
    (:g3d, :Φ),
    (:g3d, :strain_history),
    (:gh, :basal_water_thickness),
    (:gh, :hydraulic_potential_b),
    (:gh, :effective_pressure),
    (:gh, :θ_ave),
    (:gh, :preBfactor),
)

# Extra rheology fields at the end of a timestep (not part of InitialConditions).
const CHECKPOINT_DERIVED_FIELD_PATHS = (
    (:gh, :quad_f0),
    (:gh, :quad_f1),
    (:gh, :quad_f2),
    (:gh, :β),
    (:gh, :βeff),
    (:gh, :ηav),
    (:gu, :βeff),
    (:gv, :βeff),
)

const CHECKPOINT_FIELD_PATHS = (
    CHECKPOINT_INITIAL_CONDITION_FIELD_PATHS...,
    CHECKPOINT_DERIVED_FIELD_PATHS...,
)

"""
Follow `field_path` through nested fields.

For example, `(:gh, :h)` returns `fields.gh.h`.
"""
function gridfield_get(fields, field_path::Tuple)
    value = fields
    for name in field_path
        value = getproperty(value, name)
    end
    return value
end

# Build the symbol path used by MPI gather/scatter helpers, e.g. [:global_fields, :gh, :h].
checkpoint_global_path(field_path::Tuple) = Symbol[:global_fields, field_path...]

function gather_checkpoint_fields!(model::AbstractModel{<:Any, <:Any, <:MPISpec})
    for field_path in CHECKPOINT_FIELD_PATHS
        collect_mpi_field!(model, checkpoint_global_path(field_path))
    end
    return nothing
end

function copy_checkpoint_fields!(
    dst_fields,
    src_fields,
    field_paths = CHECKPOINT_FIELD_PATHS,
)
    for field_path in field_paths
        gridfield_get(dst_fields, field_path) .= gridfield_get(src_fields, field_path)
    end
    return nothing
end

"""
    basic_model_from_gathered_fields(model)

Build a serial model on the full grid from fields already gathered to root.

Copies the same physics objects as the live MPI model so the usual serial
checkpoint writer can save a portable file.
"""
function basic_model_from_gathered_fields(model::AbstractModel{<:Any, <:Any, <:MPISpec})
    gf = model.spec.global_fields
    phys = model.spec.global_physics
    phys isa NamedTuple || error(
        "MPISpec.global_physics is unset; cannot write a portable checkpoint",
    )
    serial_model = Model(
        model.spec.global_grid,
        copy(gf.gh.b),
        BasicSpec();
        initial_conditions = initial_conditions_from_fields(gf),
        params = deepcopy(phys.params),
        solver_params = deepcopy(phys.solver_params),
        shelf_melt_rate = deepcopy(phys.shelf_melt_rate),
        surface_mass_balance = deepcopy(phys.surface_mass_balance),
        fracture = deepcopy(phys.fracture),
        sliding_law = deepcopy(phys.sliding_law),
        basal_hydrology = deepcopy(phys.basal_hydrology),
        thermo_dynamics = deepcopy(phys.thermo_dynamics),
        verbose = false,
    )
    # Initial-condition fields are already set; copy derived rheology only.
    copy_checkpoint_fields!(serial_model.fields, gf, CHECKPOINT_DERIVED_FIELD_PATHS)
    # Bring wavelets and related fields in line with the restored velocities
    # so a BasicSpec pickup matches a continuous run.
    update_velocities_on_h_grid!(serial_model)
    update_shelf_strain_rate!(serial_model)
    update_surf_speed!(serial_model)
    update_model_wavelets!(serial_model)
    update_rheological_operators!(serial_model)
    return serial_model
end

function write_mpi_checkpoint!(
    model::AbstractModel{<:Any, <:Any, <:MPISpec},
    timestepping_params::TimesteppingParams,
    output_params::OutputParams,
    clock::Clock,
)
    gather_checkpoint_fields!(model)
    @root begin
        serial_model = basic_model_from_gathered_fields(model)
        write_checkpoint!(serial_model, timestepping_params, output_params, clock)
    end
    return nothing
end

"""
    finalize_checkpoint_restore!(model)

Finish an MPI pickup after fields have been scattered.

Exchange halo cells, rebuild operators, and refresh derived fields that are not
stored in the portable checkpoint but are needed for the next velocity solve.

Wavelets matter after a restart. During a normal run they are rebuilt from the
latest velocities after each velocity solve, so the next step starts with
wavelets that match those velocities. A freshly built MPI model still has the
empty wavelets from construction. If those are left in place, the multigrid
solver uses a bad starting guess. With few Picard iterations that error grows
quickly and viscosity drifts away from a continuous run. Rebuilding wavelets
here (from the restored velocities) avoids that.
"""
function finalize_checkpoint_restore!(model::AbstractModel{<:Any, <:Any, <:MPISpec})
    if model.spec.pou
        halo_exchange!(model; fields = [:h, :β, :βeff, :ηav])
    else
        halo_exchange!(model; fields = [:h, :u, :v, :β, :βeff, :ηav])
    end
    update_rheological_operators!(model)
    update_velocities_on_h_grid!(model)
    update_shelf_strain_rate!(model)
    update_surf_speed!(model)
    update_model_wavelets!(model)
    return nothing
end

"""
Load a global BasicSpec checkpoint on root, scatter its fields into the live
MPI model, finish the restore, and return a Clock on every rank.

Root copies the loaded arrays into `global_fields`. Each rank then receives the
full arrays and copies its own local patch (including halo cells).
"""
function restore_mpi_from_checkpoint!(
    model::AbstractModel{<:Any, <:Any, <:MPISpec},
    timestepping_params::TimesteppingParams,
    output_params::OutputParams,
)
    comm = model.spec.comm
    n_iter_buf = zeros(Int, 1)
    time_buf = zeros(Float64, 2)

    if model.spec.rank == 0
        loaded, clock = load_checkpoint(timestepping_params, output_params)
        copy_checkpoint_fields!(model.spec.global_fields, loaded.fields)
        n_iter_buf[1] = clock.n_iter
        time_buf[1] = clock.time
        time_buf[2] = clock.ref_time
    end
    MPI.Bcast!(n_iter_buf, comm)
    MPI.Bcast!(time_buf, comm)

    for field_path in CHECKPOINT_FIELD_PATHS
        mpi_fill_local_from_global!(
            model,
            checkpoint_global_path(field_path),
            gridfield_get(model.fields, field_path),
        )
    end
    finalize_checkpoint_restore!(model)

    return Clock(
        n_iter = n_iter_buf[1],
        time = time_buf[1],
        ref_time = time_buf[2],
    )
end

function WAVI.Simulations.pickup_model(
    model::Model{T, N, S},
    timestepping_params::TimesteppingParams,
    output_params::OutputParams,
) where {T, N, S <: MPISpec}
    clock = restore_mpi_from_checkpoint!(model, timestepping_params, output_params)
    model.spec.rank == 0 && println("Pickup successful")
    return (model, clock)
end
