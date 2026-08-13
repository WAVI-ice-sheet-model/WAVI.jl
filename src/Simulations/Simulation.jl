module Simulations

export Simulation

using JLD2
using Parameters
using Setfield
using ImageFiltering: centered, imfilter, reflect, Fill

using WAVI: AbstractModel, AbstractSimulation
using WAVI.Outputs: OutputParams, load_checkpoint
using WAVI.Parameters: TimesteppingParams
using WAVI.Time

struct Simulation{M,TS,O,C} <: AbstractSimulation
    model::M
    timestepping_params::TS
    output_params::O
    clock::C
end

"""
    Simulation(;
            model,
            timestepping_params,
            output_params = OutputParams())

Construct a WAVI.jl Simulation object.

Keyword arguments
=================
- `model`: (required) an instance of a `Model` object. When `timestepping_params.niter0 > 0`, state is
  loaded from the checkpoint file. For serial/`ThreadedSpec` runs the passed-in `model` is replaced by the
  loaded model; for `MPISpec` the driver-built model is kept and fields are scattered into it.
- `timestepping_params`: (required) an instance of a `TimesteppingParams` object, which stores information relating to timestepping
- `output_params`: an instance of an `OutputParams` object, which stores information relating to outputting of solutions.
  Also used with `timestepping_params` to locate checkpoint files for pickup (see `checkpoint_path` in `Outputs`).
"""
function Simulation(; 
                    model::AbstractModel,
                    timestepping_params::TimesteppingParams,
                    output_params::OutputParams = OutputParams())
    isnothing(timestepping_params) && throw(ArgumentError("You must specify a timestepping parameters"))

    #compute number of timesteps per output (should be robust for Inf output frequency)
    output_params = set_n_iter_out!(output_params, timestepping_params.dt, timestepping_params.n_iter_total)
    pickup_model, pickup_clock = pickup!(model, timestepping_params, output_params)

    if ~isnothing(pickup_model)
        model, clock = pickup_model, pickup_clock
    else
        # TODO: is the change from the default relevant - time is now type-variant (Real not Int)
        clock = Clock(n_iter = 0, time = 0.0, ref_time = timestepping_params.ref_time)
    end
    #set the timestep in model parameters (fudge to allow model to see the timestep in velocity solve)
    model = set_dt_in_model!(model, timestepping_params.dt)

    return Simulation(model, timestepping_params, output_params, clock)
end

Simulation(m::AbstractModel, tp::TimesteppingParams; kwargs...) = Simulation(; model=m, timestepping_params=tp, kwargs...)

include("run_simulation.jl")

function set_dt_in_model!(model, dt)
    # TODO: code smell, this should be in Model construction, requires model recreation via SetField and ConstructionBase
    model = @set model.params.dt = dt
    return model
end


function set_n_iter_out!(output_params, dt, n_iter_total)
    # TODO: code smell, this should be in the constructor for OutputParams
    output_params.output_freq == Inf ? n_iter_out = (n_iter_total + 1) : n_iter_out = round(Int, output_params.output_freq/dt)
    output_params = @set output_params.n_iter_out = n_iter_out
    return output_params
end

function pickup!(model::AbstractModel, timestepping_params::TimesteppingParams, output_params::OutputParams)
    if timestepping_params.niter0 > 0
        @info "detected niter0 > 0 (niter0 = $(timestepping_params.niter0)). Looking for pickup..."
        try
            return pickup_model(model, timestepping_params, output_params)
        catch e
            @error "Pickup error: $e"
            error("Pickup error, terminating run")
        end
    end
    return (nothing, nothing)
end

"""
    pickup_model(model, timestepping_params, output_params) -> (model, clock)

Load a checkpoint and return the model and clock to continue from.

For serial and ThreadedSpec runs, the checkpoint file replaces the model you
passed in. The returned model and clock come from that file.

For MPISpec, the fields are copied into the live distributed model instead.
The same model object is returned, together with the restored clock.
"""
function pickup_model(model::AbstractModel, timestepping_params::TimesteppingParams, output_params::OutputParams)
    loaded, clock = load_checkpoint(timestepping_params, output_params)
    println("Pickup successful")
    return (loaded, clock)
end

end
