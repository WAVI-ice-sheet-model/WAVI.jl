struct InversionData{T <: Real, N <: Integer}
    udi:: N
    vdi:: N
    hdi:: N
    speed_u::Array{T,2}
    speed_v::Array{T,2} 
    speed_u_mask::Array{Bool,2}     
    speed_v_mask::Array{Bool,2}     
    dhdt::Array{T,2}
    accumulation_rate::Array{T,2}
    dhdtaccmask::Array{Bool,2}     
end


#= """
    InversionData(; 
    nxh, nyh, nxu, nyu, nxv, nyv,
    speed_u, speed_v, dhdt, accumulation_rate)

Constructs an `InversionData` object with given data arrays.

Keyword arguments:
=================
- `nxh`, `nyh`: Dimensions of the `h` grid.
- `nxu`, `nyu`: Dimensions of the `u` grid.
- `nxv`, `nyv`: Dimensions of the `v` grid.
- `speed_u`: `u` component of surface speed (Array{T,2}).
- `speed_v`: `v` component of surface speed (Array{T,2}).
- `dhdt`: Surface height rate of change (Array{T,2}).
- `accumulation_rate`: Accumulation rate data (Array{T,2}).

Throws:
=======
- `DimensionMismatch` if any array does not match the expected grid size.
"""
function InversionData(; 
    speed_u, 
    speed_v, 
    dhdt, 
    accumulation_rate
) 

    # Check sizes of inputs
    #size(speed_u) == (nxu, nyu) || throw(DimensionMismatch("speed_u size must be (nxu x nyu), got $(size(speed_u))"))
   # size(speed_v) == (nxv, nyv) || throw(DimensionMismatch("speed_v size must be (nxv x nyv), got $(size(speed_v))"))
   # size(dhdt) == (nxh, nyh) || throw(DimensionMismatch("dhdt size must be (nxh x nyh), got $(size(dhdt))"))
   # size(accumulation_rate) == (nxh, nyh) || throw(DimensionMismatch("accumulation_rate size must be (nxh x nyh), got $(size(accumulation_rate))"))

    # Return the constructed object
    return InversionData(
                        speed_u, speed_v, dhdt, accumulation_rate)
end
 =#
function InversionData(; udi = nothing,
    vdi = nothing,
    hdi = nothing,
    speed_u = nothing, 
    speed_v = nothing,
    speed_u_mask = nothing, 
    speed_v_mask = nothing,
    dhdt = nothing,
    accumulation_rate = nothing,
    dhdtaccmask = nothing)

return InversionData(
    udi,
    vdi,
    hdi,    
    speed_u, 
    speed_v,
    speed_u_mask,
    speed_v_mask,
    dhdt,
    accumulation_rate,
    dhdtaccmask
    )
end
