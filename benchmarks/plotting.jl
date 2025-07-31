function create_heatmap_animation(netcdf_file::String, variable_name::String, output_dir::String)
    @info "Creating heatmap animation for variable: $variable_name"
    
    # Create output directory if it doesn't exist
    mkpath(output_dir)
    
    # Open NetCDF file
    nc = NCDatasets(netcdf_file)
    
    try
        # Read coordinates
        x = nc["x"]
        y = nc["y"]
        time = nc["time"]
        
        # Read the variable data
        data = nc[variable_name]
        
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
                aspect_ratio=:equal,
                size=(600, 500)
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
function plot_multiple_variables(netcdf_file::String, variables::Vector{String}, output_dir::String)
    for var in variables
        var_output_dir = joinpath(output_dir, var)
        try
            create_heatmap_animation(netcdf_file, var, var_output_dir)
        catch e
            @info "Error processing variable $var: $e"
        end
    end
end