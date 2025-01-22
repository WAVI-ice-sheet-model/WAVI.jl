struct InversionDataUGrid{T <: Real, N <: Integer}
                   nxu :: N                                    # Number of frid cells in x-direction in InversionDataUGrid
                   nyu :: N                                    # Number of grid cells in y-direction in InversionDataUGrid 
                 mask :: Array{Bool,2}                         # Mask specifying model domain wrt U grid 
           mask_inner :: Array{Bool,2}                         # Mask specifying interior of model domain wrt U grid 
                    n :: N                                     # Total number of cells in model domain 
                   ni :: N                                     # Total number of cells in interior of model domain 
                 crop :: Diagonal{T,Array{T,1}}                # Crop matrix: diagonal matrix with mask entries on diag
                 samp :: SparseMatrixCSC{T,N}                  # Sampling matrix: take full domain to model domain 
           samp_inner :: SparseMatrixCSC{T,N}                  # Sampling matrix: take full domain to interior of model domain 
               spread :: SparseMatrixCSC{T,N}                  # Spread matrix: take model domain to full domain
         spread_inner :: SparseMatrixCSC{T,N}                  # Spread matrix: take interior of model domain to full domain
              speed_u :: Array{T,2}                            # Surface speed data on u
end
    
"""
    InversionDataUGrid(;
            nxu,
            nyu,
            mask = trues(nxu,nyu),
            u_isfixed = falses(nxu,nyu)
            )

Construct a WAVI.jl InversionDataUGrid with size (nxu,nyu)
InversionDataUGrid stores fields that are defined on the problem's U grid. 
(Co-ordinates of InversionDataUGrid stored in a Grid under xxu and yyu fields)

Keyword arguments
=================
    - 'nxu': (required) Number of grid cells in x-direction in InversionDataUGrid (should be same as grid.nx + 1)
            Note that we store the grid size here, even though it can be easily inferred from grid, to increase transparency in velocity solve.
    - 'nyu': (required) Number of grid cells in y-direction in InversionDataUGrid (should be same as grid.ny)
    - 'mask': Mask specifying the model domain with respect to U grid
    - 'u_isfixed' Mask specifying where u velocities are fixed.
    - 'speed_u' surface speed on u data
"""
function InversionDataUGrid(;
                nxu,
                nyu,
                mask = trues(nxu,nyu),
                u_isfixed = falses(nxu,nyu),
                speed_u = zeros(nxu,nyu)
                )

    #check the sizes of inputs
    (size(mask) == (nxu,nyu)) || throw(DimensionMismatch("Sizes of inputs to InversionDataUGrid must all be equal to nxu x nyu (i.e. $nxu x $nyu)"))
    (size(u_isfixed) == (nxu,nyu)) || throw(DimensionMismatch("Sizes of inputs to InversionDataUGrid must all be equal to nxu x nyu (i.e. $nxu x $nyu)"))
    (size(speed_u) == (nxu,nyu)) || throw(DimensionMismatch("Sizes of inputs to InversionDataUGrid must all be equal to nxu x nyu (i.e. $nxu x $nyu)"))

    #construct operators
    n = count(mask)
    mask_inner = mask .& .! u_isfixed
    ni = count(mask_inner)
    crop = Diagonal(float(mask[:]))
    samp = sparse(1:n,(1:(nxu*nyu))[mask[:]],ones(n),n,nxu*nyu)
    samp_inner = sparse(1:ni,(1:(nxu*nyu))[mask_inner[:]],ones(ni),ni,nxu*nyu)
    spread = sparse(samp')
    spread_inner = sparse(samp_inner')


    #size assertions
    @assert n == count(mask)
    @assert ni == count(mask_inner)
    @assert crop == Diagonal(float(mask[:]))
    @assert samp == sparse(1:n,(1:(nxu*nyu))[mask[:]],ones(n),n,nxu*nyu)
    @assert samp_inner == sparse(1:ni,(1:(nxu*nyu))[mask_inner[:]],ones(ni),ni,nxu*nyu)
    @assert spread == sparse(samp')
    @assert spread_inner == sparse(samp_inner')

    #make sure boolean type rather than bitarray
    mask = convert(Array{Bool,2}, mask)
    mask_inner = convert(Array{Bool,2}, mask_inner)
    u_isfixed = convert(Array{Bool,2}, u_isfixed)

    return InversionDataUGrid(
                nxu,
                nyu,
                mask,
                mask_inner,
                n,
                ni,
                crop,
                samp,
                samp_inner,
                spread,
                spread_inner,
                speed_u
                )
end