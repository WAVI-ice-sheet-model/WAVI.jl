module Simulations

export Simulation

using JLD2
using Parameters
using Setfield

using WAVI: AbstractModel, AbstractSimulation
using WAVI.Outputs: OutputParams
using WAVI.Parameters: TimesteppingParams
using WAVI.Time

mutable struct Simulation{M,TS,O,C} <: AbstractSimulation
    model::M
    timestepping_params::TS
    output_params::O
    clock::C
end

"""
    Simulation(;
            model = nothing,
            timestepping_params = nothing,
            output_params = OutputParams(),
            pickup_output_update_flag = false)

Construct a WAVI.jl Simulation object.

Keyword arguments
=================
- `model`: (required) an instance of a `Model`` object
- `timestepping_params`: (required) an instance of a `TimesteppingParams` object, which stores information relating to timestepping
- `output_params`: an instance of an `OutputParams` object, which stores information relating to outputting of solutions
- `pickup_output_update_flag`: a flag which specifies whether to update the output_params upon picking up.
"""
function Simulation(; 
                    model::AbstractModel,
                    timestepping_params::TimesteppingParams,
                    output_params::OutputParams = OutputParams(),
                    pickup_output_update_flag = false)

    (timestepping_params !== nothing) || throw(ArgumentError("You must specify a timestepping parameters"))

    #compute number of timesteps per output (should be robust for Inf output frequency)
    output_params = set_n_iter_out!(output_params,timestepping_params.dt, timestepping_params.n_iter_total)

    #initialize the clock
    clock = Clock(n_iter = 0, time = 0.0)

    #set the timestep in model parameters (fudge to allow model to see the timestep in velocity solve)
    model = set_dt_in_model!(model, timestepping_params.dt)

    #build the simulation
    simulation = Simulation(model, timestepping_params, output_params, clock)

    #pickup?
    pickup!(simulation, pickup_output_update_flag)

    return simulation    
end

Simulation(m::AbstractModel, tp::TimesteppingParams; kwargs...) = Simulation(; model=m, timestepping_params=tp, kwargs...)

include("run_simulation.jl")

function set_dt_in_model!(model, dt)
    model = @set model.params.dt = dt
    return model
end


function set_n_iter_out!(output_params, dt,n_iter_total)
    output_params.output_freq == Inf ? n_iter_out = (n_iter_total + 1) : n_iter_out = round(Int,output_params.output_freq/dt)
    output_params = @set output_params.n_iter_out = n_iter_out
    return output_params
end

function pickup!(simulation, pickup_output_update_flag)
    #@unpack model, timestepping_params, clock = simulation

    if simulation.timestepping_params.niter0 > 0
        n_iter_string =  lpad(simulation.timestepping_params.niter0, 10, "0"); #filename as a string with 10 digits
        println("detected niter0 > 0 (niter0 = $(simulation.timestepping_params.niter0)). Looking for pickup...")
        try 
            sim_load = load(string("PChkpt_",n_iter_string, ".jld2"), "simulation")
            println("Pickup successful")

            simulation.model = sim_load.model
            simulation.clock = sim_load.clock
            simulation.output_params = sim_load.output_params #pointer based outputting system means we have to use same output parameters after pickup
        catch 
            Throw(error("Pickup error, terminating run"))
        end
    end
    return simulation
end
    
end