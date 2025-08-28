export run_simulation!, timestep!, update_clock!, update_thickness!

using WAVI
import WAVI: AbstractModel, AbstractSimulation
using WAVI.Outputs: write_outputs, zip_output
using WAVI.Processes: update_state!
using WAVI.Specs

"""
    timestep!(model, output_params, clock, timestepping_params)

Perform one timestep of the simulation
"""
function timestep!(model::AbstractModel,
                   timestepping_params::TimesteppingParams,
                   output_params::OutputParams,
                   clock::Clock)
    update_state!(model, clock)

    #write solution if at the first timestep (hack for https://github.com/RJArthern/WAVI.jl/issues/46 until synchronicity is fixed)
    # Have made the interface consistent
    # Have also removed the dependence on individual calls
    if (output_params.output_start) && (clock.n_iter == 0)
        write_outputs(model, timestepping_params, output_params, clock)
    end
    
    if timestepping_params.step_thickness
        update_thickness!(model, timestepping_params)
    end
    update_clock!(clock, timestepping_params)

    write_outputs(model, timestepping_params, output_params, clock)
end
timestep!(s::AbstractSimulation) = timestep!(s.model, s.timestepping_params, s.output_params, s.clock)

"""
update_thickness!(model::AbstractModel)

Update thickness using rate of change of thickness and apply minimum thickness constraint. Includes an option for not evolving shelf thickness.
"""
function update_thickness!(model::AbstractModel, timestepping_params)
    hUpdate = zeros(model.grid.nx,model.grid.ny)
    aground = zeros(model.grid.nx,model.grid.ny)
    hUpdate[model.fields.gh.mask] = max.(
        model.params.minimum_thickness .- model.fields.gh.h[model.fields.gh.mask],
        timestepping_params.dt * model.fields.gh.dhdt[model.fields.gh.mask])
    
    #Specify whether to evolve the shelves:
    if !model.params.evolveShelves
        hUpdate[model.fields.gh.mask] = max.(
            model.params.smallHAF .- (
                model.params.density_ocean ./ model.params.density_ice
            ) .* model.fields.gh.b[model.fields.gh.mask] .- model.fields.gh.h[model.fields.gh.mask], 
            hUpdate[model.fields.gh.mask])
        aground=(model.fields.gh.haf.>=0)
        wc=[1 1 1; 1 1 1; 1 1 1]
        w=centered(wc)
        nearfloat_mask = imfilter(model.fields.gh.mask.&.!aground,reflect(w),Fill(0,w))
        nearfloat_mask = iszero.(iszero.(nearfloat_mask))
        hUpdate[nearfloat_mask].=0
    end
    hUpdate[model.fields.gh.h_isfixed] .= 0
    model.fields.gh.h[model.fields.gh.mask] = model.fields.gh.h[model.fields.gh.mask] .+ hUpdate[model.fields.gh.mask]
end
update_thickness!(simulation::Simulation) = update_thickness!(s.model, s.timestepping_params)


"""
    update_clock!(simulation::AbstractSimulation)

Update the simulation clock
"""
function update_clock!(clock::Clock, timestepping_params::TimesteppingParams)
    clock.n_iter += 1
    clock.time += timestepping_params.dt
end
update_clock!(s::Simulation) = update_clock!(s.clock, s.timestepping_params)

"""
    run_simulation!(simulation)
    
Perform the simulation specified by the simulation
"""
function run_simulation!(simulation)
    @unpack model, timestepping_params, output_params = simulation

    for i = (simulation.clock.n_iter+1):timestepping_params.n_iter_total
        @info "Running iteration $(simulation.clock.n_iter)/$(timestepping_params.n_iter_total)"
        timestep!(simulation)
    end

    #zip the simulation output (no zipping handled by zip_output)
    zip_output(simulation)
end

