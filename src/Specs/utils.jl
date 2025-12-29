"""
ramp(i, n, overlap, leave_lower, leave_upper)

Returns a weight in [0, 1] for the partition-of-unity. The weight ramps from 0 at
the edge to 1 in the interior over `overlap` cells.

Arguments:
    i:          Index
    n:          Total number of cells
    overlap:    Number of cells to ramp over
    leave_lower(Bool): Flag. If set, do not ramp at the lower edge (i=1)
    leave_upper(Bool): Flag. If set, do not ramp at the upper edge (i=n)

Returns:
    w:          Weight in [0, 1] for cell i

Behavior:
    - If `leave_lower=false`, the weight ramps linearly from 0 → 1 over the first `overlap` cells.
    - If `leave_upper=false`, the weight ramps linearly from 1 → 0 over the last `overlap` cells.
    - If `leave_lower=true` or `leave_upper=true`, the ramp is not applied at that edge, and the weight remains 1.
"""
function ramp(i::Integer, n::Integer, overlap::Integer, leave_lower::Bool, leave_upper::Bool) :: Float64
    w = 1.0

    if !leave_lower
        w = clamp(i/overlap, 0.0, 1.0)
    end
    if !leave_upper
        w = clamp((n-i)/overlap, 0.0, 1.0)
    end

    return w
end


"""
partition_of_unity(m,n,leavei1,leaveim,leavej1,leavejn,overlapi,overlapj)

Returns a partition of unity array pou(i,j) of size m x n. The partition of unity ramps from 1 in interior to zero at edges.

    m:          Number of cells in the i-direction
    n:          Number of cells in the j-direction
    leavei1(Bool):  Flag. If set, we leave the partition of unity as one towards the edge i=1
    leaveim(Bool):  Flag. If set, we leave the partition of unity as one towards the edge i=m
    leavej1(Bool):  Flag. If set, we leave the partition of unity as one towards the edge j=1
    leavejn(Bool):  Flag. If set, we leave the partition of unity as one towards the edge j=m
    overlapi:       The ramp is applied over overlapi cells in direction i.
    overlapj:       The ramp is applied over overlapj cells in direction j.

"""
function partition_of_unity(m,n,leavei1,leaveim,leavej1,leavejn,overlapi,overlapj) :: Array{Float64,2}
    @assert overlapi ≥ 1 && overlapj ≥ 1
    @assert overlapi < m && overlapj < n
    @info "$(m) x $(n) with $(overlapi) halo in x and $(overlapj) halo in y"

    wi = [ramp(i, m, overlapi, leavei1, leaveim) for i in 1:m]
    wj = [ramp(j, n, overlapj, leavej1, leavejn) for j in 1:n]
    pou = wi * wj'

    return pou
end

"""
validate_dimension(name::String, global_dim::Int, procs::Int, halo::Int)

Checks if the local grid size is sufficient for halo overlap.

    name:         Name of the dimension (e.g., "nx", "ny")
    global_dim:   Global grid size in the dimension
    procs:        Number of processes in the dimension
    halo:         Halo size in the dimension
"""

function validate_dimension(name::String, global_dim::Int, procs::Int, halo::Int)
    procs == 1 && return

    core::Int = div(global_dim, procs)

    # If procs == 1, only have one process, req is 0
    # If procs == 2, only have edges, there is no interior
    # If procs > 2, have interior and edge nodes
    has_interior::Bool = procs > 2
    req::Int = has_interior ? 2 * halo : halo

    if core <= req
        min_global_dim = procs * (req + 1)
        max_procs = div(global_dim, req + 1)

        throw(ArgumentError(
            "MPI Configuration Error: Local grid size in $name is too small for halo overlap!\n\n" *
            "Current configuration:\n" *
            "\tGlobal $name        = $(global_dim)\n" *
            "\tProcesses        = $(procs)\n" *
            "\tHalo size        = $(halo)\n" *
            "\tLocal core size  = $(core)\n\n" *
            "Partition of Unity requires:\n" *
            "\t• Local core size > Halo (edge nodes)\n" *
            "\t• Local core size > 2*Halo (interior nodes)\n\n" *
            "Fix:\n" *
            "\t• Increase $name to at least $(min_global_dim) (keeping processes fixed), or\n" *
            "\t• Decrease processes in this direction to ≤ $(max_procs) (keeping $name fixed).\n"
        ))
    end

end
