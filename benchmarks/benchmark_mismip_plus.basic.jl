using Dates
using Printf

include("MISMIP_PLUS.jl")
include("plotting.jl")
include("utils.jl")

function main()
    @info "JULIA WAVI BENCHMARK SUITE"
    @info "Start time: $(now())"
    
    # Configuration
    output_dir = "benchmark_output_$(Dates.format(now(), "yyyymmdd_HHMMSS"))"
    mkpath(output_dir)

    # Run the model with benchmarking
    @info "Executing model..."
    result, benchmark_results = monitor_resources(MISMIP_PLUS)
    
    # Display results
    @info "Execution time: $(@sprintf("%.3f", benchmark_results.execution_time)) seconds"
    @info "Memory usage: $(@sprintf("%.2f", benchmark_results.memory_usage / 1024^2)) MB"
    @info "GC time: $(@sprintf("%.3f", benchmark_results.gc_time)) seconds"
    @info "Allocations: $(benchmark_results.allocations)"
    @info "Profile samples: $(benchmark_results.profile_data["profile_samples"])"
    
    benchmark_file = joinpath(output_dir, "benchmark_results.json")
    save_benchmark_results(benchmark_results, benchmark_file)
    
    profile_file = joinpath(output_dir, "profile_data.txt")
    open(profile_file, "w") do io
        Profile.print(io, format=:flat)
    end
    @info "Profile data saved to: $profile_file"
    
    netcdf_output = "output.nc"
    
    if isfile(netcdf_output)
        @info "Creating visualizations..."
        variables_to_plot = [:h, :u, :v, :b]
        
        plot_output_dir = joinpath(output_dir, "plots")
        plot_multiple_variables(netcdf_output, variables_to_plot, plot_output_dir)
    else
        @info "Warning: NetCDF output file '$netcdf_output' not found. Skipping visualization."
    end
    
    @info "End time: $(now())"
    @info "All output saved to: $output_dir"
end



if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
