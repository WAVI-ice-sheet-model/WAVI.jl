using NCDatasets
using Plots
using Printf

function plot_thing(ds, nom, nt; dir=".")
    lims = (minimum(ds[nom]), maximum(ds[nom]))
    for i in range(1, nt)
        filename = "$(dir)/$(string(nom)).$(@sprintf("%010.2f", ds[:TIME][i])).png"
        @info "Plotting $(nom) - t $(filename)"
        heatmap(ds[nom][:, :, i], clims=lims, c=:bluesreds)
        plot!(fmt=:png)
        savefig(filename)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    nc_file = popfirst!(ARGS)
    @info "Opening $(nc_file) with remaining args $(ARGS)"
    
    NCDataset(nc_file, "r") do ds
        @info "Data contains $(keys(ds))"
        vars = length(ARGS) > 0 ? ARGS[:] : [k for k in keys(ds) if k ∉ ("TIME", "x", "y")]
        vars = [Symbol(str) for str in vars]
        for nom in vars
            plot_thing(ds, nom, length(ds[:TIME]); dir=dirname(nc_file))
        end
    end
end