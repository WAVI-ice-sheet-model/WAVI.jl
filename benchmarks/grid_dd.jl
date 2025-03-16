using ArgParse
using BenchmarkTools
using LinearAlgebra
using Logging
using Parameters
using WAVI

@kwdef struct GridParams
    nx::Int64
    ny::Int64
    mx::Int64
    my::Int64
    cell_spacing::Float64 = 8000.0
    depth::Int64 = 4
end

function parse_cli()
    arg_settings = ArgParseSettings(exit_after_help = true)
    @add_arg_table! arg_settings begin
        "--depth", "-d"
            help = "Vertical depth on grid"
            default = 4
            arg_type = Int
        "--cell-spacing", "-s"
            help = "Override the test grid cell spacing"
            default = "8000.0"            
        "--verbose", "-v"
            help = "Turn logging debug messages on"
            action = :store_true
        "nx"
            arg_type = Int
        "ny"
            arg_type = Int
        "mx"
            default = 1
            arg_type = Int
        "my"
            default = 1
            arg_type = Int
    end
    return parse_args(ARGS, arg_settings)
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
    
    #@info "Creating model with spec"
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
    @info "Updating state for model"
    update_state!(model)
end

args = parse_cli()
@debug args
params = GridParams(
    nx = args["nx"], 
    ny = args["ny"],
    mx = args["mx"], 
    my = args["my"]
)
run_grid_ops(create_model(params, BasicParallelSpec()))
run_grid_ops(create_model(params, SharedMemorySpec()))
