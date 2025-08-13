
export write_output, write_vel

using JLD2
using MAT

using WAVI
using WAVI: AbstractModel, AbstractSimulation
using WAVI.Specs

#file containing outputting functions
"""
    write_output(model, output_params, clock)

Output the data from the simulation at the current timestep
"""
function write_output(model::AbstractModel, output_params, clock)
    output_dict = fetch_output(output_params)
    @info summary(output_dict)

    #put the grid co-ordinates and time into output.
    #Round time in output to some decimal places to make it prettier (machine precision can make this look nasty!)
    if ~haskey(output_dict, :t); output_dict["t"] = round(clock.time, digits = 3); end
    if ~haskey(output_dict, :x); output_dict["x"] = model.grid.xxh; end
    if ~haskey(output_dict, :y); output_dict["y"] = model.grid.yyh; end

    @root begin
        fname = string(output_params.output_path, output_params.prefix , lpad(clock.n_iter, 10,"0"));
        if output_params.output_format == "jld2"
            fname = string(fname, ".jld2")
            save(fname, output_dict)
        elseif output_params.output_format == "mat"
            fname = string(fname, ".mat")
            matwrite(fname, output_dict)
        end
        @info "Output at timestep number $(clock.n_iter) - $(fname)"
    end
end
write_output(s::AbstractSimulation) = write_output(s.model, s.output_params, s.clock)

"""
    function write_vel(simulation)

Write the velocity at the the final timestep of the simulation (used in the coupled wavi-mitgcm model to communicate with streamice)
"""
function write_vel(model::AbstractModel, output_params::OutputParams)
    @root begin
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
end
write_vel(s::AbstractSimulation) = write_vel(s.model, s.output_params)


