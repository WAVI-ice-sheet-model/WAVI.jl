"""
partition_of_unity(m,n,leavei1,leaveim,leavej1,leavejn,overlapi,overlapj)

Returns a partition of unity array pou(i,j) of size m x n. The partition of unity ramps from 1 in interior to zero at edges.

leavei1(Bool):  Flag. If set, we leave the partition of unity as one towards the edge i=1
leaveim(Bool):  Flag. If set, we leave the partition of unity as one towards the edge i=m
leavej1(Bool):  Flag. If set, we leave the partition of unity as one towards the edge j=1 
leavejn(Bool):  Flag. If set, we leave the partition of unity as one towards the edge j=m  
overlapi:       The ramp is applied over overlapi cells in direction i.
overlapj:       The ramp is applied over overlapj cells in direction j.

"""
function partition_of_unity(m,n,leavei1,leaveim,leavej1,leavejn,overlapi,overlapj)
    @assert m > (~leavei1 && overlapi) + (~leaveim && overlapi)
    @assert n > (~leavej1 && overlapj) + (~leavejn && overlapj)
    pou = [min(1.0, 
        leavei1 ? Inf : (i-1)./(overlapi),
        leaveim ? Inf : (m-i)./(overlapi)) for i=1:m, j=1:n] .*
        [min(1.0, 
        leavej1 ? Inf : (j-1)./(overlapj),
        leavejn ? Inf : (n-j)./(overlapj)) for i=1:m, j=1:n]
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
