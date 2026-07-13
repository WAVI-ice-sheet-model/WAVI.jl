using CSV
using NCDatasets
using Plots
using Printf

"""Return a NetCDF variable, matching `name` exactly or case-insensitively."""
function nc_var(ds, name::AbstractString)
    haskey(ds, name) && return ds[name]
    lname = lowercase(name)
    for k in keys(ds)
        lowercase(k) == lname && return ds[k]
    end
    throw(KeyError(name))
end

function create_heatmap_animation(netcdf_file::String, variable_name::String, output_dir::String)
    @info "Creating heatmap animation for variable: $variable_name"
    
    # Create output directory if it doesn't exist
    mkpath(output_dir)
    
    # Open NetCDF file
    nc = NCDataset(netcdf_file)
    
    try
        # Read coordinates (WAVI output may use mixed-case names, e.g. TIME)
        x = nc_var(nc, "x")
        y = nc_var(nc, "y")
        time = nc_var(nc, "time")
        data = nc_var(nc, variable_name)

        @info "Data dimensions: $(size(data))"
        @info "Time steps: $(length(time))"
        
        # Create heatmaps for each time step
        for (i, t) in enumerate(time)
            # Extract 2D slice for this time step
            slice_2d = data[:, :, i]
            
            # Create heatmap
            p = heatmap(
                x, y, slice_2d',
                title="$variable_name at time = $(@sprintf("%.2f", t))",
                xlabel="X",
                ylabel="Y",
                color=:viridis,
                clim=(minimum(data), maximum(data)),
                aspect_ratio=:equal,
                size=(600, 500),
            )
            
            # Save frame
            frame_filename = joinpath(output_dir, "$(variable_name)_frame_$(lpad(i, 4, '0')).png")
            savefig(p, frame_filename)
            
            if i % 10 == 0 || i == length(time)
                @info "Processed frame $i/$(length(time))"
            end
        end
        
        @info "Animation frames saved to: $output_dir"
        
    finally
        close(nc)
    end
end

"""
    plot_multiple_variables(netcdf_file::String, variables::Vector{String}, output_dir::String)

Create heatmap animations for multiple variables from NetCDF file.
"""
function plot_multiple_variables(netcdf_file::String, 
                                 variables::Vector{String}, 
                                 output_dir::String)
    for var in variables
        var_output_dir = joinpath(output_dir, var)
        try
            create_heatmap_animation(netcdf_file, var, var_output_dir)
        catch e
            @info "Error processing variable $var: $e"
        end
    end
end

"""
    plot_resource_timeseries(csv_paths; labels=String[], reference_cores=nothing, output=...)

Overlay RSS and CPU time series from one or more `resource_timeseries.csv` files.

CPU is plotted as a fraction: uses a `cpu_fraction` column if present, otherwise
`cpu_cores_used / reference_cores`. When `reference_cores` is omitted, each CSV
picks it up from the run’s `benchmark_results.json` (`metadata.reference_cores`),
falling back to `1`. Saves a two-subplot PNG.
"""
function plot_resource_timeseries(
    csv_paths::Vector{<:AbstractString};
    labels::Vector{<:AbstractString} = String[],
    reference_cores::Union{Nothing, Real} = nothing,
    output::AbstractString = joinpath(BENCHMARK_OUTPUT_DIR, "resource_comparison.png"),
)
    isempty(csv_paths) && error("Pass at least one resource_timeseries.csv path.")
    reference_cores === nothing || reference_cores > 0 || error("reference_cores must be positive.")

    labs = if isempty(labels)
        [basename(dirname(abspath(p))) for p in csv_paths]
    else
        collect(String, labels)
    end
    length(labs) == length(csv_paths) ||
        error("Number of labels ($(length(labs))) must match number of CSV paths ($(length(csv_paths))).")

    plt_rss = Plots.plot(; xlabel = "elapsed (s)", ylabel = "RSS (MB)", legend = :topright)
    plt_cpu = Plots.plot(; xlabel = "elapsed (s)", ylabel = "CPU fraction", legend = :topright)

    for (path, lab) in zip(csv_paths, labs)
        isfile(path) || error("CSV not found: $path")
        rows = CSV.File(path)
        t = collect(rows.elapsed_s)
        rss_mb = collect(rows.rss_bytes) ./ 1024^2
        cols = propertynames(rows)
        frac = if :cpu_fraction in cols
            collect(rows.cpu_fraction)
        else
            ref = something(reference_cores, _reference_cores_for_csv(path))
            collect(rows.cpu_cores_used) ./ Float64(ref)
        end
        Plots.plot!(plt_rss, t, rss_mb; label = lab)
        Plots.plot!(plt_cpu, t, frac; label = lab)
    end

    fig = Plots.plot(plt_rss, plt_cpu; layout = (2, 1), size = (800, 600))
    out_dir = dirname(output)
    !isempty(out_dir) && mkpath(out_dir)
    Plots.savefig(fig, output)
    @info "Saved resource comparison to: $output"
    return output
end

"""Read `metadata.reference_cores` from the run’s `benchmark_results.json`, or `1`."""
function _reference_cores_for_csv(csv_path::AbstractString)
    json_path = joinpath(dirname(abspath(csv_path)), "benchmark_results.json")
    isfile(json_path) || return 1.0
    try
        data = JSON3.read(read(json_path, String))
        meta = get(data, :metadata, nothing)
        if meta !== nothing && haskey(meta, :reference_cores)
            return Float64(meta.reference_cores)
        end
        haskey(data, :reference_cores) && return Float64(data.reference_cores)
    catch e
        @debug "Could not read reference_cores from $json_path" exception = e
    end
    return 1.0
end
