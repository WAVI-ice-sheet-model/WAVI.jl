mutable struct InversionSimulation{M1,M2,TS,O,C} 
    model::M1
    inversion:: M2
    JKVstepping_params::TS
    output_params::O
    clock::C
end

"""
    InversionSimulation(;
            model = nothing,
            JKVstepping_params = nothing,
            output_params = OutputParams(),
            pickup_output_update_flag = false)

Construct a WAVI.jl InversionSimulation object.

Keyword arguments
=================
- `model`: (required) an instance of a `Model`` object
- `inversion`: (required) an instance of a `Model`` object
- `JKVstepping_params`: (required) an instance of a `TimesteppingParams` object, which stores information relating to JKVstepping
- `output_params`: an instance of an `OutputParams` object, which stores information relating to outputting of solutions
- `pickup_output_update_flag`: a flag which specifies whether to update the output_params upon picking up.
"""
function InversionSimulation(;
                    model = nothing,
                    inversion = nothing,
                    JKVstepping_params = nothing,
                    output_params = OutputParams(),
                    pickup_output_update_flag = false)

    (JKVstepping_params !== nothing) || throw(ArgumentError("You must specify a JKVstepping parameters"))

   # #set number of JKVsteps per output:
   # output_params.n_iter_out =  JKVstepping_params.n_iter_out

    #initialize the clock
    clock = Clock(n_iter = 0, time = 0.0)

   # #set the JKVstep in model parameters (fudge to allow model to see the JKVstep in velocity solve)
   #  model = set_dt_in_model!(model, JKVstepping_params.dt)

    #build the simulation
    inversion_simulation = InversionSimulation(model, inversion, JKVstepping_params, output_params, clock)

    #pickup? NOT SETUP FOR PICKUP YET!!!
    pickup_inversion!(inversion_simulation, pickup_output_update_flag)

    return inversion_simulation 
    
end

include("run_inversion_simulation.jl")

 function pickup_inversion!(inversion_simulation, pickup_output_update_flag)
    #@unpack model, JKVstepping_params, clock = simulation

    if inversion_simulation.JKVstepping_params.niter0 > 0
        n_iter_string =  lpad(inversion_simulation.JKVstepping_params.niter0, 10, "0"); #filename as a string with 10 digits
        println("detected niter0 > 0 (niter0 = $(inversion_simulation.JKVstepping_params.niter0)). Looking for pickup...")
        try 
            #for running on vscode vis workstation:
            #filename = string(inversion_simulation.output_params.output_path, "/PChkpt_", n_iter_string, ".jld2")
            #for running with ensembler:
            filename = string("PChkpt_", n_iter_string, ".jld2")
            println("Looking for file: ", filename)

            if isfile(filename)
                println("file exists")
                sim_load = load(filename, "inversion_simulation")
                println("Pickup successful")
            else
                error("File $filename not found!")
            end
            
            inversion_simulation.model = sim_load.model
            inversion_simulation.inversion = sim_load.inversion
            inversion_simulation.clock = sim_load.clock
            inversion_simulation.output_params = sim_load.output_params #pointer based outputting system means we have to use same output parameters after pickup

        catch e
            println("Error while loading: ", e)
            rethrow()  # Re-throws the original error so you can see it
        end
    end
    return inversion_simulation
end
    