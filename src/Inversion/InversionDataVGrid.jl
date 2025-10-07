struct InversionDataVGrid{T <: Real, N <: Integer}
                   nxv :: N                                    # Number of frid cells in x-direction in InversionDataVGrid
                   nyv :: N                                    # Number of grid cells in y-direction in InversionDataVGrid 
                 mask :: Array{Bool,2}                         # Mask specifying model domain wrt U grid 
           mask_inner :: Array{Bool,2}                         # Mask specifying interior of model domain wrt U grid 
                    n :: N                                     # Total number of cells in model domain 
                   ni :: N                                     # Total number of cells in interior of model domain 
                 crop :: Diagonal{T,Array{T,1}}                # Crop matrix: diagonal matrix with mask entries on diag
                 samp :: SparseMatrixCSC{T,N}                  # Sampling matrix: take full domain to model domain 
           samp_inner :: SparseMatrixCSC{T,N}                  # Sampling matrix: take full domain to interior of model domain 
               spread :: SparseMatrixCSC{T,N}                  # Spread matrix: take model domain to full domain
         spread_inner :: SparseMatrixCSC{T,N}                  # Spread matrix: take interior of model domain to full domain
              speed_v :: Array{T,2}                            # Surface speed data on v
end
    
"""
    InversionDataVGrid(;
                nxv,
                nyv,
                mask = trues(nxv,nyv),
                v_isfixed = falses(nxv,nyv),
                speed_v = zeros(nxv,nyv)
                )

Construct a WAVI.jl InversionDataVGrid with size (nxv,nyv)
InversionDataVGrid stores fields that are defined on the problem's v grid. 
(Co-ordinates of InversionDataVGrid stored in a Grid under xxu and yyu fields)

Keyword arguments
=================
    - 'nxv': (required) Number of grid cells in x-direction in InversionDataVGrid (should be same as grid.nx + 1)
            Note that we store the grid size here, even though it can be easily inferred from grid, to increase transparency in velocity solve.
    - 'nyv': (required) Number of grid cells in y-direction in InversionDataVGrid (should be same as grid.ny)
    - 'mask': Mask specifying the model domain with respect to v grid
    - 'u_isfixed': Mask specifying where v velocities are fixed.
    - 'speed_v': surface speed on v data
"""
function InversionDataVGrid(;
                nxv,
                nyv,
                mask = trues(nxv,nyv),
                v_isfixed = falses(nxv,nyv),
                speed_v = zeros(nxv,nyv)
                )

#    @unpack gu,gv,gh = model.fields

    #check the sizes of inputs
    (size(mask) == (nxv,nyv)) || throw(DimensionMismatch("Sizes of inputs to InversionDataVGrid must all be equal to nxv x nyv (i.e. $nxv x $nyv)"))
    (size(v_isfixed) == (nxv,nyv)) || throw(DimensionMismatch("Sizes of inputs to InversionDataVGrid must all be equal to nxv x nyv (i.e. $nxv x $nyv)"))
    (size(speed_v) == (nxv,nyv)) || throw(DimensionMismatch("Sizes of inputs to InversionDataVGrid must all be equal to nxv x nyv (i.e. $nxv x $nyv)"))


    #refine masks to only include points in the model mask:
   # mask = ((mask .== 1) .& (gv.mask .== 1))
    mask=convert(Array{Bool,2}, mask)


    #construct operators
    n = count(mask)
    mask_inner = mask .& .! v_isfixed
    ni = count(mask_inner)
    crop = Diagonal(float(mask[:]))
    samp = sparse(1:n,(1:(nxv*nyv))[mask[:]],ones(n),n,nxv*nyv)
    samp_inner = sparse(1:ni,(1:(nxv*nyv))[mask_inner[:]],ones(ni),ni,nxv*nyv)
    spread = sparse(samp')
    spread_inner = sparse(samp_inner')

#=     #Do smoothing here:
    vs_data_vec=speed_v[mask]
    vs_data_spread=spread*vs_data_vec
    #
    vs_data_crop = gv.crop * vs_data_spread[:]  
  #  vs_data_crop = gv.crop * vdata[:]  
    vs_data_cent = gv.cent * vs_data_crop 
    vs_data_cent_samp = gh.samp  * vs_data_cent
    vs_data_cent_samp_spread= gh.spread*vs_data_cent_samp
   #
    vs_data_centT =  gv.centᵀ*vs_data_cent_samp_spread
    vs_data_centT_crop =  gv.crop*vs_data_centT
    #
    vs_data_sampi = samp_inner*vs_data_centT_crop

    vs_data_smoothed=zeros(gv.nxv,gv.nyv)
    vs_data_smoothed[mask]=vs_data_sampi
    speed_v=vs_data_smoothed

    ##pre-select points not near data gaps:  
 #   v_neargap = imfilter((gv.mask .& .!(mask .> 0)) .|> Int, Kernel.ones(3,3), Pad(1,1)) .> 0
     data_gap = gv.mask .& .!mask 
     outside_model = .!gv.mask
     gap_mask = data_gap .| outside_model 
      v_neargap = imfilter((gap_mask) .|> Int, Kernel.ones(3,3),  Pad(:replicate)) .> 0
    mask = ((mask .== 1) .& (gv.mask .== 1) .& (v_neargap .== 0)) 


     #  mask[1, :] .= false
     #  mask[end, :] .= false
     #  mask[:,1] .= false
      # mask[:, end] .= false

    #construct operators
    n = count(mask)
    mask_inner = mask .& .! v_isfixed
    ni = count(mask_inner)
    crop = Diagonal(float(mask[:]))
    samp = sparse(1:n,(1:(nxv*nyv))[mask[:]],ones(n),n,nxv*nyv)
    samp_inner = sparse(1:ni,(1:(nxv*nyv))[mask_inner[:]],ones(ni),ni,nxv*nyv)
    spread = sparse(samp')
    spread_inner = sparse(samp_inner') =#


    #size assertions
    @assert n == count(mask)
    @assert ni == count(mask_inner)
    @assert crop == Diagonal(float(mask[:]))
    @assert samp == sparse(1:n,(1:(nxv*nyv))[mask[:]],ones(n),n,nxv*nyv)
    @assert samp_inner == sparse(1:ni,(1:(nxv*nyv))[mask_inner[:]],ones(ni),ni,nxv*nyv)
    @assert spread == sparse(samp')
    @assert spread_inner == sparse(samp_inner')

    #make sure boolean type rather than bitarray
    mask = convert(Array{Bool,2}, mask)
    mask_inner = convert(Array{Bool,2}, mask_inner)
    v_isfixed = convert(Array{Bool,2}, v_isfixed)

    return InversionDataVGrid(
                nxv,
                nyv,
                mask,
                mask_inner,
                n,
                ni,
                crop,
                samp,
                samp_inner,
                spread,
                spread_inner,
                speed_v
                )
end