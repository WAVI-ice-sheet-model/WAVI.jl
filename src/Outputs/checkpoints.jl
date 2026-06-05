using JLD2

using WAVI.Parameters: TimesteppingParams
using WAVI.Time: Clock

"""
    checkpoint_path(timestepping_params, output_params)

Directory for permanent checkpoint files (`Chkpt_*.jld2`).

Uses `timestepping_params.chkpt_path` when set explicitly; otherwise falls back to
`output_params.output_path` so runs that only set an output folder still write and
read checkpoints in the same place.
"""
function checkpoint_path(timestepping_params::TimesteppingParams, output_params::OutputParams)
    if timestepping_params.chkpt_path != "./"
        return timestepping_params.chkpt_path
    elseif output_params.output_path != "./"
        return output_params.output_path
    else
        return timestepping_params.chkpt_path
    end
end

checkpoint_filename(n_iter::Integer) = string("Chkpt_", lpad(n_iter, 10, "0"), ".jld2")

function should_write_checkpoint(timestepping_params::TimesteppingParams, clock::Clock)
    return timestepping_params.chkpt_freq != Inf &&
           clock.n_iter > 0 &&
           mod(clock.n_iter, timestepping_params.n_iter_chkpt) == 0
end

function write_checkpoint!(model, timestepping_params::TimesteppingParams, output_params::OutputParams, clock::Clock)
    path = checkpoint_path(timestepping_params, output_params)
    if !isdir(path)
        mkpath(path)
    end
    fname = joinpath(path, checkpoint_filename(clock.n_iter))
    @save fname model=model timestepping_params=timestepping_params clock=clock
    @info "Permanent checkpoint at timestep number $(clock.n_iter) — $(fname)"
    return nothing
end

function load_checkpoint(timestepping_params::TimesteppingParams, output_params::OutputParams)
    path = checkpoint_path(timestepping_params, output_params)
    n_iter = timestepping_params.niter0
    fname = joinpath(path, checkpoint_filename(n_iter))
    dict = load(fname)
    if haskey(dict, "simulation")
        # legacy format (pre #105): @save fname simulation saved everything
        # in one bundle, need to unpack it.
        sim_load = dict["simulation"]
        model = sim_load.model
        clock = sim_load.clock
    else
        # new format (post #105): @save fname model=model clock=clock saved separately
        model = dict["model"]
        clock = dict["clock"]
    end
    if clock.n_iter != n_iter
        @warn "Checkpoint iteration $(clock.n_iter) does not match niter0=$n_iter ($(fname))"
    end
    return model, clock
end
