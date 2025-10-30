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