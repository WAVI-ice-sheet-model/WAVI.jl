using WAVI

struct MPISpec
    top::Integer
    right::Integer
    bottom::Integer
    left::Integer
    halo::Integer
    global_grid_nx::Integer
    global_grid_ny::Integer
    px::Integer
    py::Integer
    coords::Array{Integer, 1}
    global_size::Integer
    rank::Integer
end

# Emulate MPI_PROC_NULL (typically -1 in MPI.jl)
const PROC_NULL = -1

function emulate_mpi_creation(global_size::Integer, px::Integer, py::Integer, current_rank::Integer)
    """
    # Testing

    # # Usage: Emulate rank 0 in a 2x1 grid (global_size=2)
    # # (x=0, y=0) -> neighbors should be Right=1, others PROC_NULL
    coords, top, right, bottom, left = emulate_mpi_creation(2, 2, 1, 0)

    # Usage: Emulate a grid with size 6 (e.g., 2x3 grid) and get neighbors for rank 1
    # coords, top, right, bottom, left = emulate_mpi_creation(6, 3, 2, 1)
    """

    # Calculate coordinates (x, y) for the given rank
    # Note: MPI Cart_create with dims=(px, py) varies the last dimension (y) fastest
    x_coord = div(current_rank, py)  # Column number
    y_coord = current_rank % py      # Row number
    coords = [x_coord, y_coord]

    # Function to calculate the destination rank based on shift direction
    function get_neighbor(direction, displacement)
        # Apply shift
        nx, ny = x_coord, y_coord
        if direction == 0  # X-direction (left/right)
            nx += displacement
        else  # Y-direction (up/down)
            ny += displacement
        end

        # Check boundaries and return destination rank (or PROC_NULL if out of bounds)
        if nx < 0 || nx >= px || ny < 0 || ny >= py
            return PROC_NULL
        else
            return nx * py + ny  # Rank calculation based on column-major order (x*py + y)
        end
    end

    # Emulate Cart_shift for all 4 directions (top, right, bottom, left)
    top    = get_neighbor(1, -1)  # Y-direction, shift -1 (Up)
    right  = get_neighbor(0, 1)   # X-direction, shift +1 (Right)
    bottom = get_neighbor(1, 1)   # Y-direction, shift +1 (Down)
    left   = get_neighbor(0, -1)  # X-direction, shift -1 (Left)

    return coords, top, right, bottom, left
end


"""
Get the halo size for each direction

Returns (top, right, bottom, left) halo sizes

    Args:
        halo (Integer): Halo size
    Returns:
        Tuple{Int, Int, Int, Int}: Halo size for each direction

"""
function get_halos(spec::MPISpec)::Tuple{Int, Int, Int, Int}
    # Note: MPI_Cart_shift returns `MPI_PROC_NULL` (represented as a negative integer)
    #       if the neighbour is not present
    return spec.top > -1 ? spec.halo : 0,
           spec.right > -1 ? spec.halo : 0,
           spec.bottom > -1 ? spec.halo : 0,
           spec.left > -1 ? spec.halo : 0
end

# Get size of local grid
function get_size(spec::MPISpec)::Tuple{Int, Int}
    nx_base = div(spec.global_grid_nx, spec.px)
    ny_base = div(spec.global_grid_ny, spec.py)
    nx_rem  = rem(spec.global_grid_nx, spec.px)
    ny_rem  = rem(spec.global_grid_ny, spec.py)

    # If this rank's x-coord is < remainder, it gets an extra cell
    # This accounts for the fact that the global grid may not be evenly split
    # (e.g. if the global grid is 10x10 and we have 3 ranks in x, then the first
    # rank will have 4 cells in x and the other two ranks will have 3 cells in x)
    local_nx = nx_base + (spec.coords[1] < nx_rem ? 1 : 0)
    local_ny = ny_base + (spec.coords[2] < ny_rem ? 1 : 0)

    halos = get_halos(spec)
    local_nx += halos[4] + halos[2]
    local_ny += halos[1] + halos[3]

    return local_nx, local_ny
end

# Get bounds of local grid
function get_bounds(spec::MPISpec)::Tuple{Int, Int, Int, Int}
    # Calculate base local grid size assuming even distribution
    nx_base = div(spec.global_grid_nx, spec.px)
    nx_rem  = rem(spec.global_grid_nx, spec.px)
    ny_base = div(spec.global_grid_ny, spec.py)
    ny_rem  = rem(spec.global_grid_ny, spec.py)

    # Adding extra cell to local grid size if there is a remainder
    local_nx = nx_base + (spec.coords[1] < nx_rem ? 1 : 0)
    local_ny = ny_base + (spec.coords[2] < ny_rem ? 1 : 0)

    # Calculate start index (accounting for previous ranks that took an extra cell)
    # (Rank * BaseSize) + (How many previous ranks took an extra cell) + 1
    x_start = spec.coords[1] * nx_base + min(spec.coords[1], nx_rem) + 1
    y_start = spec.coords[2] * ny_base + min(spec.coords[2], ny_rem) + 1
    
    # Calculate end index
    x_end = x_start + local_nx - 1
    y_end = y_start + local_ny - 1

    halos = get_halos(spec)
    return max(x_start - halos[4], 1),                   # Left: Protect against indices < 1
           min(x_end   + halos[2], spec.global_grid_nx), # Right: Protect against indices > global max
           max(y_start - halos[1], 1),                   # Bottom: Protect against indices < 1
           min(y_end   + halos[3], spec.global_grid_ny)  # Top: Protect against indices > global max
end

function MISMIP_PLUS_GRID(;
        nx = 80,
        ny = 10,
    )
    #Grid and boundary conditions
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
    return grid
end
