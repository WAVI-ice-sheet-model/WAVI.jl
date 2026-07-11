# # Resource monitor — memory and CPU samples for benchmarks
#
# On Linux this starts a second, small Julia process. That helper reads the
# target process's memory and CPU use from `/proc/<pid>/` every few tenths of a
# second and writes the results to a CSV file.
#
# On other OS, only record peak memory when sampling stops (`Sys.maxrss()`).
#
# Command to run self-test:
#   julia --project=benchmarks -t 1 benchmarks/resource_monitor.jl
#
# Helper (started by `start_monitor!`; do not run by hand unless debugging):
#   julia --startup-file=no benchmarks/resource_monitor.jl --sampler \
#       <pid> <interval_s> <csv_path> <stop_path>

module ResourceMonitorModule

export ResourceMonitor, start_monitor!, stop_monitor!
export read_rss, read_peak_rss_hwm, normalise_cpu_fraction

"""
    ResourceSample

One row of resource telemetry.

# Fields
- `elapsed_s` — seconds since `start_monitor!` began
- `rss_bytes` — how much RAM the process is using right now (Linux: `VmRSS`)
- `cpu_cores_used` — CPU effort over the last sample gap, as a number of cores
  (about `1.0` means one core was fully busy; more than `1` means several
  threads or BLAS workers were busy at once)
"""
struct ResourceSample
    elapsed_s::Float64
    rss_bytes::Int
    cpu_cores_used::Float64
end

"""
    ResourceMonitor

Handle returned by `start_monitor!`. Holds the sample list, peak memory, and
the paths/process used by the external helper. Prefer `start_monitor!` /
`stop_monitor!` over building this yourself.
"""
mutable struct ResourceMonitor
    samples::Vector{ResourceSample}
    peak_rss_bytes::Int
    sample_interval_s::Float64
    active::Bool
    process::Union{Nothing, Base.Process}
    csv_path::String
    stop_path::String
    tmpdir::String
    start_ns::UInt64
end

# How many CPU time units equal one second on this machine (Linux `CLK_TCK`).
const CLK_TCK = Sys.islinux() ? ccall(:sysconf, Clong, (Cint,), 2) : 100

const _SCRIPT_PATH = @__FILE__

"""
    read_rss(pid = getpid()) -> Int

Current resident memory (RSS) for process `pid`, in bytes.

On Linux this reads `VmRSS` from `/proc/<pid>/status`. Elsewhere it falls back
to `Sys.maxrss()` (peak so far for this Julia process, not a live reading).
"""
function read_rss(pid::Integer = getpid())
    Sys.islinux() || return Int(Sys.maxrss())
    try
        for line in eachline("/proc/$pid/status")
            if startswith(line, "VmRSS:")
                m = match(r"(\d+)\s+kB", line)
                m !== nothing && return parse(Int, m.captures[1]) * 1024
            end
        end
    catch e
        @debug "Failed to read /proc/$pid/status" exception = e
    end
    return 0
end

"""
    read_peak_rss_hwm(pid = getpid()) -> Int

Highest resident memory seen so far for process `pid`, in bytes.

On Linux this reads `VmHWM` (“high water mark”) from `/proc/<pid>/status`.
Elsewhere it falls back to `Sys.maxrss()`.
"""
function read_peak_rss_hwm(pid::Integer = getpid())
    Sys.islinux() || return Int(Sys.maxrss())
    try
        for line in eachline("/proc/$pid/status")
            if startswith(line, "VmHWM:")
                m = match(r"(\d+)\s+kB", line)
                m !== nothing && return parse(Int, m.captures[1]) * 1024
            end
        end
    catch e
        @debug "Failed to read VmHWM from /proc/$pid/status" exception = e
    end
    return 0
end

"""
    read_cpu_ticks(pid = getpid()) -> Int

Total CPU time used by process `pid`, in kernel clock ticks.

This is user time plus system time (`utime` + `stime` in `/proc/<pid>/stat`).
Convert to seconds by dividing by `CLK_TCK`. Returns `0` off Linux or on error.
"""
function read_cpu_ticks(pid::Integer = getpid())
    Sys.islinux() || return 0
    try
        stat_str = read("/proc/$pid/stat", String)
        # The process name sits in parentheses and may contain spaces, so find
        # the last ')' and only then split the numeric fields that follow.
        idx = findlast(')', stat_str)
        idx === nothing && return 0
        parts = split(stat_str[idx + 2:end])
        # After ') ': field 12 = utime, field 13 = stime (1-based in this slice).
        return parse(Int, parts[12]) + parse(Int, parts[13])
    catch e
        @debug "Failed to read /proc/$pid/stat" exception = e
        return 0
    end
end

"""
    normalise_cpu_fraction(cpu_cores_used, reference_cores) -> Float64

Turn a “cores used” reading into a fraction of the cores you meant to use.

For example, if the run was allowed 4 threads and `cpu_cores_used` is `2.0`,
the fraction is `0.5`. Pass Julia threads, Schwarz subdomain count, or SLURM
tasks as `reference_cores`, depending on the mode. Returns `0` if
`reference_cores` is not positive.
"""
function normalise_cpu_fraction(cpu_cores_used::Real, reference_cores::Real)
    reference_cores > 0 || return 0.0
    return Float64(cpu_cores_used) / Float64(reference_cores)
end

"""
    run_external_sampler!(pid, interval_s, csv_path, stop_path)

Helper entry point: sample process `pid` every `interval_s` seconds and append
rows to `csv_path`.

Stops when `stop_path` is created, or when `/proc/<pid>` disappears (the target
process has exited). Intended to run in the child process started by
`start_monitor!`.
"""
function run_external_sampler!(pid::Int, interval_s::Float64, csv_path::String, stop_path::String)
    start_ns = time_ns()
    prev_ticks = read_cpu_ticks(pid)
    prev_ns = start_ns

    open(csv_path, "w") do io
        println(io, "elapsed_s,rss_bytes,cpu_cores_used")
        rss0 = read_rss(pid)
        println(io, join(("0.0", string(rss0), "0.0"), ","))
        flush(io)

        while !isfile(stop_path)
            sleep(interval_s)
            isfile(stop_path) && break
            isdir("/proc/$pid") || break

            now_ns = time_ns()
            rss = read_rss(pid)
            ticks = read_cpu_ticks(pid)
            dt = (now_ns - prev_ns) / 1e9
            cpu_cores = dt > 0 ? ((ticks - prev_ticks) / CLK_TCK) / dt : 0.0
            elapsed = (now_ns - start_ns) / 1e9
            println(io, join((string(elapsed), string(rss), string(cpu_cores)), ","))
            flush(io)

            prev_ticks = ticks
            prev_ns = now_ns
        end
    end
    return nothing
end

"""Load `ResourceSample` rows from helper CSV (skips header line)."""
function _load_samples_csv(csv_path::String)
    samples = ResourceSample[]
    isfile(csv_path) || return samples
    open(csv_path, "r") do io
        readline(io)  # header
        for line in eachline(io)
            isempty(strip(line)) && continue
            parts = split(line, ',')
            length(parts) >= 3 || continue
            push!(
                samples,
                ResourceSample(
                    parse(Float64, parts[1]),
                    parse(Int, parts[2]),
                    parse(Float64, parts[3]),
                ),
            )
        end
    end
    return samples
end

"""
    start_monitor!(interval_s = 0.25) -> ResourceMonitor

Begin recording memory and CPU use for **this** process.

On Linux, starts the external helper and waits until it has written the first
sample (so Julia start-up time does not leave a blank gap at the beginning).
Sampling does not depend on `Threads.nthreads()`.

On other platforms, returns a monitor that only fills in peak memory later in
`stop_monitor!` (no time series).

Call `stop_monitor!` when the work you want to measure has finished.
"""
function start_monitor!(interval_s::Float64 = 0.25)
    mon = ResourceMonitor(ResourceSample[], 0, interval_s, true, nothing, "", "", "", time_ns())
    Sys.islinux() || return mon

    tmpdir = mktempdir(; prefix = "wavi_rm_")
    csv_path = joinpath(tmpdir, "samples.csv")
    stop_path = joinpath(tmpdir, "stop")
    log_path = joinpath(tmpdir, "sampler.log")

    pid = getpid()
    cmd = `$(Base.julia_cmd()) --startup-file=no --history-file=no $_SCRIPT_PATH --sampler $pid $interval_s $csv_path $stop_path`
    try
        proc = run(pipeline(cmd; stdout = log_path, stderr = log_path); wait = false)
        mon.process = proc
        mon.csv_path = csv_path
        mon.stop_path = stop_path
        mon.tmpdir = tmpdir

        # Wait until the helper has started and written the t = 0 sample.
        deadline = time() + 60.0
        while time() < deadline
            if !process_running(proc)
                log_txt = isfile(log_path) ? read(log_path, String) : ""
                @warn "ResourceMonitor: sampler exited before first sample" log = log_txt
                break
            end
            if isfile(csv_path) && count(==('\n'), read(csv_path, String)) >= 2
                break
            end
            sleep(0.05)
        end
    catch e
        @warn "ResourceMonitor: failed to start external sampler; peak RSS only" exception = e
        try
            rm(tmpdir; recursive = true, force = true)
        catch
        end
        mon.active = false
        mon.process = nothing
        mon.csv_path = ""
        mon.stop_path = ""
        mon.tmpdir = ""
    end
    return mon
end

"""
    stop_monitor!(mon) -> ResourceMonitor

Stop sampling and finalise results on `mon`.

Writes the stop file so the helper exits, waits for it, loads the CSV into
`mon.samples`, then sets `mon.peak_rss_bytes` to the best peak available
(largest sample, Linux `VmHWM`, and current RSS — or `Sys.maxrss()` off Linux).
Removes the temporary directory afterwards.
"""
function stop_monitor!(mon::ResourceMonitor)
    mon.active = false

    if !isempty(mon.stop_path)
        try
            open(mon.stop_path, "w") do io
                write(io, "stop\n")
            end
        catch e
            @debug "Failed to write sampler stop file" exception = e
        end
    end

    if mon.process !== nothing
        try
            wait(mon.process)
        catch e
            @debug "Resource monitor process error on shutdown" exception = e
        end
        mon.process = nothing
    end

    if !isempty(mon.csv_path)
        mon.samples = _load_samples_csv(mon.csv_path)
        if !isempty(mon.samples)
            mon.peak_rss_bytes = maximum(s.rss_bytes for s in mon.samples)
        end
    end

    if Sys.islinux()
        mon.peak_rss_bytes = max(mon.peak_rss_bytes, read_peak_rss_hwm(), read_rss())
    else
        mon.peak_rss_bytes = Int(Sys.maxrss())
    end

    if !isempty(mon.tmpdir)
        try
            rm(mon.tmpdir; recursive = true, force = true)
        catch e
            @debug "Failed to remove sampler temporary directory" exception = e
        end
        mon.tmpdir = ""
    end
    return mon
end

end # module ResourceMonitorModule

using .ResourceMonitorModule

# Child process entry: `start_monitor!` launches this file with `--sampler`.
if abspath(PROGRAM_FILE) == ResourceMonitorModule._SCRIPT_PATH && !isempty(ARGS) && ARGS[1] == "--sampler"
    length(ARGS) == 5 || error("usage: resource_monitor.jl --sampler <pid> <interval_s> <csv_path> <stop_path>")
    ResourceMonitorModule.run_external_sampler!(
        parse(Int, ARGS[2]),
        parse(Float64, ARGS[3]),
        ARGS[4],
        ARGS[5],
    )
    exit(0)
end

# Self-test when this file is run as the main script.
if abspath(PROGRAM_FILE) == ResourceMonitorModule._SCRIPT_PATH
    let
        println("Testing ResourceMonitor on $(Sys.islinux() ? "Linux" : "non-Linux")...")
        println("Julia threads: $(Threads.nthreads())")
        mon = start_monitor!(0.1)

        # Busy work that does not call `sleep` — checks that sampling still works
        # when the only Julia thread is fully occupied.
        a = zeros(Float64, 800, 800)
        t_end = time() + 2.0
        while time() < t_end
            a .+= 1.0e-6
            s = sum(a)
            a[1] = s * 1.0e-20
        end

        stop_monitor!(mon)

        println("Peak RSS: $(round(mon.peak_rss_bytes / 1024^2; digits = 2)) MB")
        println("Samples collected: $(length(mon.samples))")
        println("Sample interval: $(mon.sample_interval_s) s")
        if !isempty(mon.samples)
            max_cpu = maximum(s.cpu_cores_used for s in mon.samples)
            println("Max CPU cores used: $(round(max_cpu; digits = 2))")
            println("CPU fraction (ref=1 core): $(round(normalise_cpu_fraction(max_cpu, 1); digits = 2))")
            println("\nFirst 3 samples:")
            for (i, s) in enumerate(first(mon.samples, 3))
                println(
                    "  [$i] t=$(round(s.elapsed_s; digits = 2))s  ",
                    "RSS=$(round(s.rss_bytes / 1024^2; digits = 2))MB  ",
                    "CPU=$(round(s.cpu_cores_used; digits = 2)) cores",
                )
            end
            # About 2 s of work at 0.1 s spacing should give many samples on -t 1.
            if Sys.islinux() && length(mon.samples) < 8
                error("Expected ≥8 samples during 2s busy work; got $(length(mon.samples)) (sampler starved?)")
            end
        elseif Sys.islinux()
            error("No samples collected on Linux")
        end
        println("Test completed.")
    end
end
