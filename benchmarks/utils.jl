using Dates
using JSON3
using Profile

struct BenchmarkResults
    execution_time::Float64
    memory_usage::Float64
    gc_time::Float64
    allocations::Int64
    #profile_data::Dict{String, Any}
    system_info::Dict{String, Any}
    timestamp::DateTime
end

function benchmark_main(id::String,
                        model::Function, 
                        model_args::Dict, 
                        variables_to_plot::Vector{String},
                        rank::Int = 0)

    output_dir = "benchmark_$(id)_$(Dates.format(now(), "yyyymmdd_HHMMSS"))"
    model_args[:folder] = output_dir

    if rank == 0
        @info "JULIA WAVI BENCHMARK SUITE - $(id)"
        @info "Start time: $(now())"

        # Configuration
        mkpath(output_dir)
        @info "Executing model with $(model_args)"
    end
    
    # Run the model with benchmarking
    result, benchmark_results = monitor_resources(model; model_args...)
    
    if rank == 0
        # Display results
        @info "Execution time: $(@sprintf("%.3f", benchmark_results.execution_time)) seconds"
        @info "Memory usage: $(@sprintf("%.2f", benchmark_results.memory_usage / 1024^2)) MB"
        @info "GC time: $(@sprintf("%.3f", benchmark_results.gc_time)) seconds"
        @info "Allocations: $(benchmark_results.allocations)"
        #@info "Profile samples: $(benchmark_results.profile_data["profile_samples"])"
            
        benchmark_file = joinpath(output_dir, "benchmark_results.json")
        save_benchmark_results(benchmark_results, benchmark_file)
        
        #profile_file = joinpath(output_dir, "profile_data.txt")
        #open(profile_file, "w") do io
        #    Profile.print(io, format=:flat)
        #end
        #@info "Profile data saved to: $profile_file"
        
        netcdf_output = "$(output_dir)/outfile.nc"
        
        if isfile(netcdf_output)
            @info "Creating visualizations..."
            plot_output_dir = joinpath(output_dir, "plots")
            plot_multiple_variables(netcdf_output, variables_to_plot, plot_output_dir)
        else
            @info "Warning: NetCDF output file '$netcdf_output' not found. Skipping visualization."
        end
        
        @info "End time: $(now())"
        @info "All output saved to: $output_dir"
    end
end

function monitor_resources(func, args...; kwargs...)
    Profile.clear()
    
    @profile begin
        result = @timed func(args...; kwargs...)
    end
    
    #profile_data = Profile.fetch()
    
    benchmark_result = BenchmarkResults(
        result.time,
        result.bytes,
        result.gctime,
        result.gcstats.allocd,
        #Dict("profile_samples" => length(profile_data)),
        get_system_info(),
        now()
    )
    
    return result.value, benchmark_result
end

function get_system_info()
    return Dict(
        "julia_version" => string(VERSION),
        "cpu_threads" => Sys.CPU_THREADS,
        "total_memory" => Sys.total_memory(),
        "hostname" => gethostname(),
        "os" => string(Sys.KERNEL)
    )
end

function save_benchmark_results(results::BenchmarkResults, filename::String)
    output_data = Dict(
        "timestamp" => string(results.timestamp),
        "execution_time_seconds" => results.execution_time,
        "memory_usage_bytes" => results.memory_usage,
        "gc_time_seconds" => results.gc_time,
        "allocations" => results.allocations,
        #"profile_samples" => results.profile_data["profile_samples"],
        "system_info" => results.system_info
    )
    
    open(filename, "w") do io
        JSON3.pretty(io, output_data)
    end
    
    @info "Benchmark results saved to: $filename"
end