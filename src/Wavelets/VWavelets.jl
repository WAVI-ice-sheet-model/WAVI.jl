export VWavelets

#Struct to hold information on wavelet-grid (v-component).
struct VWavelets{T <: Real, N <: Integer}
        nxvw :: N                                       # Number of grid points in x in UWavelets (equal to UGrid)
        nyvw :: N                                       # Number of grid points in y in UWavelets (equal to UGrid)
        mask :: Array{Bool,2}                           # Model domain on the U grid
           n :: Base.RefValue{N}                        # Number of grid points in domain
        crop :: Base.RefValue{Diagonal{T,Array{T,1}}}   # Crop matrix: diagonal matrix with mask entries on diag
        samp :: Base.RefValue{SparseMatrixCSC{T,N}}     # Sampling matrix: take full domain to model domain 
      spread :: Base.RefValue{SparseMatrixCSC{T,N}}     # Spread matrix: take model domain to full domain
      levels :: N                                       # Number of wavelet levels 
        idwt :: KronType{T,N}                           # Wavelet matrix cross produce
       idwtᵀ :: KronType{T,N}                           # Adjoint of wavelet matrix cross produce
    wavelets :: Array{T,2}                              # Wavelet matrix
    correction_coarse::Base.RefValue{Vector{T} }        # cache to store coarse correction for multigrid                           
  end


"""
    VWavelets(;
            nxvw,
            nyvw,
            mask = trues(nxvw,nyvw), 
            levels)

Construct a WAVI.jl VWavelets structure with size (nxvw,nyvw)
UWavelets stores wavelet matrix and associated fields for velocities in y direction

Keyword arguments
=================
    - 'nxvw': (required) Number of grid cells in x-direction in VWavelets 
    - 'nyvw': (required) Number of grid cells in y-direction in VWavelets 
    - 'mask': Mask specifying the model domain with respect to UWavelet grid (VGrid)
    - levels: (required) Number of levels in the preconditioner 
"""
function VWavelets(;   
                    nxvw,
                    nyvw,
                    mask = trues(nxvw,nyvw),
                    levels)
     
    (size(mask) == (nxvw,nyvw)) || throw(DimensionMismatch("Sizes of inputs to UWavelets must all be equal to nxvw x nyvw (i.e. $nxvw x $nyvw)"))

    #compute non-inputs 
    n = Ref(count(mask));  @assert n[] == count(mask)
    crop = Ref(Diagonal(float(mask[:]))); @assert crop[] == Diagonal(float(mask[:]));
    samp  = Ref(sparse(1:n[],(1:(nxvw*nyvw))[mask[:]],ones(n[]),n[],nxvw*nyvw)); @assert samp[] == sparse(1:n[],(1:(nxvw*nyvw))[mask[:]],ones(n[]),n[],nxvw*nyvw)
    spread = Ref(sparse(samp[]')); @assert spread[] == sparse(samp[]')
    idwt =  wavelet_matrix(nyvw,levels,"reverse" ) ⊗ wavelet_matrix(nxvw,levels,"reverse")
    idwtᵀ =  sparse(wavelet_matrix(nyvw,levels,"reverse" )') ⊗ sparse(wavelet_matrix(nxvw,levels,"reverse")')
    wavelets = zeros(nxvw,nyvw); @assert size(wavelets)==(nxvw,nyvw)
    correction_coarse = Ref(zeros(n[])); @assert length(correction_coarse[])==n[]

    #make sure boolean type rather than bitarray
    mask = convert(Array{Bool,2}, mask)

    return VWavelets(
                    nxvw,
                    nyvw,
                    mask,
                    n,
                    crop, 
                    samp, 
                    spread, 
                    levels,
                    idwt,
                    idwtᵀ,
                    wavelets,
                    correction_coarse)
end
    