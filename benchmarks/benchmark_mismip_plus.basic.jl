using Dates
using Printf

include("MISMIP_PLUS.jl")
include("plotting.jl")
include("utils.jl")

if abspath(PROGRAM_FILE) == @__FILE__
    benchmark_main("basic", MISMIP_PLUS, Dict(), ["h", "u", "v"])
end
