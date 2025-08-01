using WAVI 

include("MISMIP_PLUS.jl")
include("plotting.jl")
include("utils.jl")

if abspath(PROGRAM_FILE) == @__FILE__
    benchmark_main("threaded", MISMIP_PLUS, Dict{Symbol, Any}(
        :spec => ThreadedSpec(ngridsx=2, ngridsy=2, overlap=2, niterations=2,)
    ), ["h", "u", "v"])
end
