using MPI
using WAVI 

include("MISMIP_PLUS.jl")
include("plotting.jl")
include("utils.jl")

if abspath(PROGRAM_FILE) == @__FILE__
    # TODO: ewww, had to clone due to requirement for presence of grid at spec creation time
    nx = 80
    ny = 10
    nσ = 4
    x0 = 0.0
    y0 = -40000.0
    dx = 8000.0
    dy = 8000.0
    h_mask=trues(nx,ny)
    u_iszero = falses(nx+1,ny); u_iszero[1,:].=true
    v_iszero=falses(nx,ny+1); v_iszero[:,1].=true; v_iszero[:,end].=true
    grid = Grid(nx = nx, 
                ny = ny,   
                nσ = nσ, 
                x0 = x0, 
                y0 = y0, 
                dx = dx, 
                dy = dy,
                h_mask = h_mask, 
                u_iszero = u_iszero, 
                v_iszero = v_iszero)
    mpi_spec = MPISpec(2, 1, 2, grid)

    benchmark_main("mpi", MISMIP_PLUS, Dict{Symbol, Any}(
        :spec => mpi_spec
    ), ["h", "u", "v"], mpi_spec.rank)
end