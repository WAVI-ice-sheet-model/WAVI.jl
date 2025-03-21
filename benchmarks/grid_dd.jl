using BenchmarkTools
using LinearAlgebra
using Logging
using Parameters
using WAVI

@kwdef struct GridParams{T, N}
    nx::T
    ny::T
    mx::T
    my::T
    cell_spacing::N = 8000.0
    depth::T = 4
end

function GridParams(nx, ny; mx=1, my=1, depth=4, cell_spacing=8000.0)
    return GridParams{Int64, Float64}(nx, ny, mx, my, depth,cell_spacing)
end

function create_model(p::GridParams, spec::AbstractParallelSpec)::WAVI.AbstractModel
    @info "Creating $(p.nx)x$(p.ny)x$(p.depth) grid"
    grid = Grid(
        nx=p.nx,
        ny=p.ny,
        dx=p.cell_spacing,
        dy=p.cell_spacing,
        nσ=p.depth,
        x0 = 0.0,
        y0 = -40000.0,
    )
    
    @info "Creating model with spec $(typeof(spec))"
    model = Model(
        grid = grid, 
        bed_elevation = [ℯ/ (y/p.ny) + 50 *sin(x * 1e-3 * p.nx) for x in range(1., 10., length = p.nx), y in range(1., 10., length = p.ny)],
        params = Params(),
        solver_params = SolverParams(),
        initial_conditions = InitialConditions(),
        melt_rate = UniformMeltRate(),
        parallel_spec = spec)

    return model
end

function run_grid_ops(model::WAVI.AbstractModel)
    @info "Updating state for model with spec $(typeof(spec))"
    update_state!(model)
end

function run_benchmark()
    run_grid_ops(create_model(params, BasicSpec()))
    run_grid_ops(create_model(params, SharedMemorySpec()))
end
