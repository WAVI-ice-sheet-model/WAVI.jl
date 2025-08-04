
export fetch_output, write_output

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
    @root begin
        output_dict = fetch_output(output_params.outputs)

        #put the grid co-ordinates and time into output.
        #Round time in output to some decimal places to make it prettier (machine precision can make this look nasty!)
        if ~haskey(output_dict, :t); output_dict["t"] = round(clock.time, digits = 3); end
        if ~haskey(output_dict, :x); output_dict["x"] = model.global_grid.xxh; end
        if ~haskey(output_dict, :y); output_dict["y"] = model.global_grid.yyh; end

        fname = string(output_params.output_path, output_params.prefix , lpad(clock.n_iter, 10,"0"));
        if output_params.output_format == "jld2"
            fname = string(fname, ".jld2")
            save(fname, output_dict)
        elseif output_params.output_format == "mat"
            fname = string(fname, ".mat")
            matwrite(fname, output_dict)
        end

        println("outputting at timestep number $(clock.n_iter)")
    end
end
write_output(s::AbstractSimulation) = write_output(s.model, s.output_params, s.clock)

"""
    fetch_output(outputs)

Return a dictionary with dictionary entries corresponding to outputs
"""
function fetch_output(outputs)
    output_dict = Dict()
    for (k,v) in zip(keys(outputs), outputs)
        output_dict[string(k)] = v
    end
    return output_dict
end

