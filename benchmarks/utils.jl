using Dates
using JSON3
using Profile

struct BenchmarkResults
    execution_time::Float64
    memory_usage::Float64
    gc_time::Float64
    allocations::Int64
    profile_data::Dict{String, Any}
    system_info::Dict{String, Any}
    timestamp::DateTime
end

function monitor_resources(func, args...; kwargs...)
    Profile.clear()
    
    @profile begin
        result = @timed func(args...; kwargs...)
    end
    
    profile_data = Profile.fetch()
    
    benchmark_result = BenchmarkResults(
        result.time,
        result.bytes,
        result.gctime,
        result.gcstats.allocd,
        Dict("profile_samples" => length(profile_data)),
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
        "profile_samples" => results.profile_data["profile_samples"],
        "system_info" => results.system_info
    )
    
    open(filename, "w") do io
        JSON3.pretty(io, output_data)
    end
    
    @info "Benchmark results saved to: $filename"
end