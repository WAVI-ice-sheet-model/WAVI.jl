
export write_outputs, write_output, write_vel

using JLD2
using MAT

using WAVI
using WAVI: AbstractModel, AbstractSimulation
using WAVI.Deferred
using WAVI.Specs


function write_outputs(model::AbstractModel, 
                       timestepping_params::TimesteppingParams, 
                       output_params::OutputParams, 
                       clock::Clock)
    #check if we have hit a permanent checkpoint
    if mod(i, timestepping_params.n_iter_chkpt) == 0
        #output a permanent checkpoint
        n_iter_string =  lpad(clock.n_iter, 10, "0"); #filename as a string with 10 digits
        fname = joinpath(output_params.output_path, string("Chkpt_",n_iter_string, ".jld2"))
        # TODO: strip out OutputParams
        @save fname simulation
        println("making permanent checkpoint at timestep number $(simulation.clock.n_iter)")
    end

    #check if we have hit an output timestep
    if mod(i,simulation.output_params.n_iter_out) == 0
        collect!(output_params, model)
        write_output(output_params, lpad(clock.n_iter, 10,"0"))
        clear!(output_params)
    end

    #check the dump velocity flag at the final timestep
    if (i == timestepping_params.n_iter_total) && output_params.dump_vel
        write_vel(output_params, model)
    end
end
write_outputs(sim::AbstractSimulation) = write_outputs(sim.model, sim.timestepping_params, sim.output_params, sim.clock)

#file containing outputting functions
"""
    write_output(model, output_params, clock)

Output the data from the simulation at the current timestep
"""
function write_output(output_params::OutputParams, name::String)
    output_dict = collect!(output_params)

    #put the grid co-ordinates and time into output.
    #Round time in output to some decimal places to make it prettier (machine precision can make this look nasty!)
    if ~haskey(output_dict, :t); output_dict["t"] = round(clock.time, digits = 3); end
    if ~haskey(output_dict, :x); output_dict["x"] = model.grid.xxh; end
    if ~haskey(output_dict, :y); output_dict["y"] = model.grid.yyh; end

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

"""
    function write_vel(simulation)

Write the velocity at the the final timestep of the simulation (used in the coupled wavi-mitgcm model to communicate with streamice)

    TODO: This can be made more generic for binary dumps of any field
"""
function write_vel(output_params::OutputParams, model::AbstractModel)
    uVel_file_string = string(output_params.prefix,  "_U.bin")
    vVel_file_string = string(output_params.prefix,  "_V.bin")
    
    u_out=model.fields.gu.u[1:end-1,:]
    v_out=model.fields.gv.v[:,1:end-1]

    u_out .= hton.(u_out)
    v_out .= hton.(v_out)

    ufileID =  open(uVel_file_string,"w")
    write(ufileID, u_out[:,:])
    close(ufileID) 
    vfileID =  open(vVel_file_string,"w")
    write(vfileID, v_out[:,:])
    close(vfileID)   
end


