struct InversionDataHGrid{T <: Real, N <: Integer}
                   nxh :: N                                    # Number of frid cells in x-direction in InversionDataHGrid
                   nyh :: N                                    # Number of grid cells in y-direction in InversionDataHGrid 
                 mask :: Array{Bool,2}                         # Mask specifying model domain wrt H grid 
            h_isfixed :: Array{Bool,2}                         # Mask specifying locations of fixed thickness
                    n :: N                                     # Total number of cells in model domain 
                 crop :: Diagonal{T,Array{T,1}}                # Crop matrix: diagonal matrix with mask entries on diag
                 samp :: SparseMatrixCSC{T,N}                  # Sampling matrix: take full domain to model domain 
               spread :: SparseMatrixCSC{T,N}                  # Spread matrix: take model domain to full domain
                  dhdt :: Array{T,2}                            # dhdt data on h
     accumulation_rate :: Array{T,2}                            # accumulation rate data on h
             residual :: Array{T,2}                            # residuals on h
end
    
"""
    InversionDataHGrid(;
            nxh,
            nyh,
            mask = trues(nxh,nyh),
             h_isfixed = falses(nxh,nxy),
            )

Construct a WAVI.jl InversionDataHGrid with size (nxh,nyh)
InversionDataHGrid stores fields that are defined on the problem's H grid. 
(Co-ordinates of InversionDataHGrid stored in a Grid under xxu and yyu fields)

Keyword arguments
=================
    - 'nxh': (required) Number of grid cells in x-direction in InversionDataHGrid (should be same as grid.nx + 1)
            Note that we store the grid size here, even though it can be easily inferred from grid, to increase transparency in velocity solve.
    - 'nyh': (required) Number of grid cells in y-direction in InversionDataHGrid (should be same as grid.ny)
    - 'mask': Mask specifying the model domain with respect to H grid
    - 'h_isfixed' Mask specifying where h velocities are fixed.
    - 'speed_u' surface speed on h data
"""
function InversionDataHGrid(;
                nxh,
                nyh,
                mask = trues(nxh,nyh),
                h_isfixed = falses(nxh,nxy),
                dhdt = zeros(nxh,nyh),
                accumulation_rate = zeros(nxh,nyh),
                residual = zeros(nxh,nyh)
                )

    #check the sizes of inputs
    (size(mask) == (nxh,nyh)) || throw(DimensionMismatch("Sizes of inputs to InversionDataHGrid must all be equal to nxh x nyh (i.e. $nxh x $nyh)"))
    (size(dhdt) == (nxh,nyh)) || throw(DimensionMismatch("Sizes of inputs to InversionDataHGrid must all be equal to nxh x nyh (i.e. $nxh x $nyh)"))
    (size(accumulation_rate) == (nxh,nyh)) || throw(DimensionMismatch("Sizes of inputs to InversionDataHGrid must all be equal to nxh x nyh (i.e. $nxh x $nyh)"))

    #construct operators
    n = count(mask)
    crop = Diagonal(float(mask[:]))
    samp = sparse(1:n,(1:(nxh*nyh))[mask[:]],ones(n),n,nxh*nyh)
    spread = sparse(samp')


    #size assertions
    @assert n == count(mask)
    @assert crop == Diagonal(float(mask[:]))
    @assert samp == sparse(1:n,(1:(nxh*nyh))[mask[:]],ones(n),n,nxh*nyh)
    @assert spread == sparse(samp')


    #make sure boolean type rather than bitarray
    mask = convert(Array{Bool,2}, mask)
    h_isfixed = convert(Array{Bool,2}, h_isfixed)

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
                residual
                )
end