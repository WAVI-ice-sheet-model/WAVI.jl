struct InversionDataHGrid{T <: Real, N <: Integer}
                   nxh :: N                                    # Number of grid cells in x-direction in InversionDataHGrid
                   nyh :: N                                    # Number of grid cells in y-direction in InversionDataHGrid 
                 mask :: Array{Bool,2}                         # Mask specifying model domain wrt H grid 
            h_isfixed :: Array{Bool,2}                         # Mask specifying locations of fixed thickness
                    n :: N                                     # Total number of cells in model domain 
                 crop :: Diagonal{T,Array{T,1}}                # Crop matrix: diagonal matrix with mask entries on diag
                 samp :: SparseMatrixCSC{T,N}                  # Sampling matrix: take full domain to model domain 
               spread :: SparseMatrixCSC{T,N}                  # Spread matrix: take model domain to full domain
                  dhdt :: Array{T,2}                            # dhdt data on h
    accumulation_rate  :: Array{T,2}                            # accumulation rate data on h
                    us :: Array{T,2}                            # u speed data (surf)  on h
                    vs :: Array{T,2}                            # v speed data (surf)  on h
            surf_speed :: Array{T,2}                            # speed data (surf)  on h
       surf_speed_mask :: Array{Bool,2}                         # surf speed data mask of valid points on h
             residual  :: Array{T,2}                            # residuals on h
end
    
"""
    InversionDataHGrid(;
                nxh,
                nyh,
                mask = trues(nxh,nyh),
                h_isfixed = falses(nxh,nyh),
                dhdt = nothing,
                accumulation_rate = zeros(nxh,nyh)
            #    us=zeros(nxh,nyh),
            #    vs=zeros(nxh,nyh),
             #   surf_speed=zeros(nxh,nyh),
             #   surf_speed_mask=falses(nxh,nyh),
             #   residual = zeros(nxh,nyh)
                )

Construct a WAVI.jl InversionDataHGrid with size (nxh,nyh)
InversionDataHGrid stores fields that are defined on the problem's data H grid. 
(Co-ordinates of InversionDataHGrid stored in a Grid under xxh and yyh fields)

Keyword arguments
=================
    - 'nxh': (required) Number of grid cells in x-direction in DataHGrid (should be same as grid.nx)
    - 'nyh': (required) Number of grid cells in y-direction in DataHGrid (should be same as grid.ny)
    - 'mask': Mask specifying the model domain
    - 'h_isfixed': Mask specifying points where ice thickness is fixed
    - 'dhdt': dhdt data 
    - 'accumulation_rate': accumlate rate data 
#    - 'us': u surface speed data
#    - 'vs': v surface speed data 
#    - 'surf_speed': surface speed data
#    - 'surf_speed_mask': mask of valid speed data
 #   - 'residual':     
"""

function InversionDataHGrid(;
                nxh,
                nyh,
                mask = trues(nxh,nyh),
                h_isfixed = falses(nxh,nyh),
                dhdt = nothing,
                accumulation_rate = zeros(nxh,nyh)
              #  us=zeros(nxh,nyh),
              #  vs=zeros(nxh,nyh),
               # surf_speed=zeros(nxh,nyh),
               # surf_speed_mask=falses(nxh,nyh),
               # residual = zeros(nxh,nyh)
                )

    #check the sizes of inputs
   
    if dhdt !== nothing
    (size(dhdt) == (nxh,nyh)) || throw(DimensionMismatch("Sizes of inputs to InversionDataHGrid must all be equal to nxh x nyh (i.e. $nxh x $nyh)"))
    else
        dhdt = dhdt === nothing ? fill(NaN, nxh, nyh) : dhdt
        println("WARNING: dhdt is not provided and is set as NaNs")
    end
    if mask !== nothing
        (size(mask) == (nxh,nyh)) || throw(DimensionMismatch("Sizes of inputs to InversionDataHGrid must all be equal to nxh x nyh (i.e. $nxh x $nyh)"))
    else 
            mask = mask === nothing ? fill(false, nxh, nyh) : mask
            println("WARNING: dhdtacc_mask is not provided: only velocities will be used in inversion")
        end
        if accumulation_rate !== nothing
    (size(accumulation_rate) == (nxh,nyh)) || throw(DimensionMismatch("Sizes of inputs to InversionDataHGrid must all be equal to nxh x nyh (i.e. $nxh x $nyh)"))
        else
                accumulation_rate = accumulation_rate === nothing ? fill(NaN, nxh, nyh) : accumulation_rate
                println("WARNING: accumulation is not provided as is set as NaNs")
        end


    #construct operators
    n = count(mask)
    crop = Diagonal(float(mask[:]))
    samp = sparse(1:n,(1:(nxh*nyh))[mask[:]],ones(n),n,nxh*nyh)
    spread = sparse(samp')

    #construct quantities not passed
    us = zeros(nxh,nyh) 
    vs = zeros(nxh,nyh)
    surf_speed = zeros(nxh,nyh) 
    surf_speed_mask = falses(nxh,nyh)
    residual = zeros(nxh,nyh)

    #size assertions
    @assert size(mask)==(nxh,nyh); #@assert mask == clip(mask)
    @assert size(h_isfixed)==(nxh,nyh); 
    @assert n == count(mask)
    @assert crop == Diagonal(float(mask[:]))
    @assert samp == sparse(1:n,(1:(nxh*nyh))[mask[:]],ones(n),n,nxh*nyh)
    @assert spread == sparse(samp')    
    @assert size(us)==(nxh,nyh)
    @assert size(vs)==(nxh,nyh)
    @assert size(surf_speed)==(nxh,nyh) 
    @assert size(surf_speed_mask)==(nxh,nyh)
    @assert size(residual)==(nxh,nyh)


    #make sure boolean type rather than bitarray
    mask = convert(Array{Bool,2}, mask)
    h_isfixed = convert(Array{Bool,2}, h_isfixed)
    surf_speed_mask = convert(Array{Bool,2}, surf_speed_mask)

    return InversionDataHGrid(
                nxh,
                nyh,
                mask,
                h_isfixed,
                n,
                crop,
                samp,
                spread,
                dhdt,
                accumulation_rate,
                us,
                vs,
                surf_speed,
                surf_speed_mask,
                residual
                )
end