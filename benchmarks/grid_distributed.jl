using BenchmarkTools
using Distributed
using LinearAlgebra
using Logging
using Parameters
using WAVI

include("grid_dd.jl");

function run_distributed()
    @info "Starting distributed.jl job"
    @warn "This implementation does not work in any way, please do not use"
    p = GridParams(100, 100; mx=2, my=2)
    dist_spec = DistributedSpec(
        ngridsx = p.mx,
        ngridsy = p.my,
        niterations = 10
    )
    model = create_model(p, dist_spec)
    run_grid_ops(model)
end

run_distributed()
