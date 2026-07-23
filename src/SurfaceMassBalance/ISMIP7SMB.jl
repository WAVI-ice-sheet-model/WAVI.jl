export ISMIP7SMB

using WAVI: AbstractClimateForcing
using NCDatasets

struct ISMIP7SMB{T <: Real,
                P <: String,
                RE <: Union{Array{T,2}, Nothing}, 
                VSG <: Union{Array{T,2}, Nothing},
                SA <: Union{Array{T,2}, Nothing},
                RS <: Union{Array{T,2}, Nothing}} <: AbstractSurfaceMassBalance
    smb_anomaly_prefix::P
    vertical_smb_gradient_prefix::P
    reference_elevation::RE
    vertical_smb_gradient::VSG
    smb_anomaly::SA
    reference_smb::RS
    path_to_forcing::String
end

function ISMIP7SMB(; 
                smb_anomaly_prefix = nothing,
                vertical_smb_gradient_prefix = nothing,
                reference_elevation = nothing,
                vertical_smb_gradient = nothing, 
                smb_anomaly = nothing,           
                reference_smb = nothing,
                path_to_forcing = "./")

    #check that you've passed a prefix for the smb anomaly and vertical smb gradient          
    ~(smb_anomaly_prefix === nothing) || throw(ArgumentError("You must pass a prefix for the smb anomaly"))
    ~(vertical_smb_gradient_prefix === nothing) || throw(ArgumentError("You must pass a prefix for the vertical smb gradient"))

    # check that you've passed reference elevation 
    ~(reference_elevation === nothing) || throw(ArgumentError("You must pass a reference elevation file"))

    #check that you've passed a vertical smb gradient
    ~(reference_smb === nothing) || throw(ArgumentError("You must pass a reference smb"))

    return ISMIP7SMB(smb_anomaly_prefix,vertical_smb_gradient_prefix,reference_elevation,vertical_smb_gradient, smb_anomaly,reference_smb, path_to_forcing)

end

function reconstruct_on_grid(smb::ISMIP7SMB,grid::Grid) 
    return ISMIP7SMB(
    smb.smb_anomaly_prefix,
    smb.vertical_smb_gradient_prefix,
    isnothing(smb.reference_elevation) ? zeros(grid.nx,grid.ny) : 
    size(smb.reference_elevation) == (grid.nx,grid.ny) ? smb.reference_elevation :
    throw(DimensionMismatch("Size of reference elevation is incompatible with grid")),
    isnothing(smb.vertical_smb_gradient) ? zeros(grid.nx,grid.ny) :
    size(smb.vertical_smb_gradient) == (grid.nx,grid.ny) ? smb.vertical_smb_gradient :
    throw(DimensionMismatch("Size of vertical smb gradient is incompatible with grid")),
    isnothing(smb.smb_anomaly) ? zeros(grid.nx,grid.ny) :
    size(smb.smb_anomaly) == (grid.nx,grid.ny) ? smb.smb_anomaly :
    throw(DimensionMismatch("Size of smb anomaly is incompatible with grid")),
    isnothing(smb.reference_smb) ? zeros(grid.nx,grid.ny) : 
    size(smb.reference_smb) == (grid.nx,grid.ny) ? smb.reference_smb :
    throw(DimensionMismatch("Size of reference smb is incompatible with grid")),
    smb.path_to_forcing)
end

function reconstruct_on_subdomain(smb::ISMIP7SMB,grid::Grid,subdomain::NTuple{4,<: Integer}) 
    x_start,x_end,y_start,y_end = subdomain
    return ISMIP7SMB(
    smb.smb_anomaly_prefix,
    smb.vertical_smb_gradient_prefix,
    size(smb.reference_elevation) == size(grid)[1:2] ? smb.reference_elevation[x_start:x_end, y_start:y_end] : smb.reference_elevation,
    size(smb.vertical_smb_gradient) == size(grid)[1:2] ? smb.vertical_smb_gradient[x_start:x_end, y_start:y_end] : smb.vertical_smb_gradient,
    size(smb.smb_anomaly) == size(grid)[1:2] ? smb.smb_anomaly[x_start:x_end, y_start:y_end] : smb.smb_anomaly,
    size(smb.reference_smb) == size(grid)[1:2] ? smb.reference_smb[x_start:x_end, y_start:y_end] : smb.reference_smb,
    smb.path_to_forcing)
end

function update_accumulation_rate!(surface_mass_balance::ISMIP7SMB, model::AbstractModel, clock::Clock)
    @unpack reference_smb, smb_anomaly, vertical_smb_gradient, reference_elevation = surface_mass_balance
    @unpack accumulation,s = model.fields.gh
    
    change_in_elevation = s - reference_elevation

    #check the files are the right size
    (size(smb_anomaly) == (model.grid.nx, model.grid.ny)) || throw(DimensionMismatch("Size of read in smb anomaly is not compatible with grid size"))
    (size(vertical_smb_gradient) == (model.grid.nx, model.grid.ny)) || throw(DimensionMismatch("Size of read in vertical smb gradient is not compatible with grid size"))


    #set the updated smb 
    accumulation .=  reference_smb .+ smb_anomaly .+ vertical_smb_gradient .* change_in_elevation
    return nothing 


end


function update_climate_forcing!(surface_mass_balance::ISMIP7SMB, grid::Grid, clock::Clock) 
    @unpack smb_anomaly = surface_mass_balance
    @unpack vertical_smb_gradient = surface_mass_balance
    @unpack dx = grid
    @unpack smb_anomaly_prefix, vertical_smb_gradient_prefix,path_to_forcing = surface_mass_balance

    #get the year from clock for the forcing files
    current_time = clock.time + clock.ref_time
    println(current_time)
    current_time_string = string(Int(round(current_time)))

    # load in the smb anomaly from ISMIP7
    resolution = join([string(Int(dx)), "m"])
    smb_anomaly_filename = joinpath(path_to_forcing,join([smb_anomaly_prefix,  current_time_string,".nc"]))
    smb_anomaly_ncfile   = NCDataset(smb_anomaly_filename)
    smb_anomaly .= replace(smb_anomaly_ncfile["acabf-anomaly"][:,:,1] , missing => NaN)
    #println("read in smb anomaly forcing file: " * smb_anomaly_filename)
    @info "read in smb anomaly forcing file: $smb_anomaly_filename"


    # load in the vertical smb gradient from ISMIP7
    vertical_smb_gradient_anomaly_filename = joinpath(path_to_forcing, join([vertical_smb_gradient_prefix,  current_time_string,".nc"]))
    vertical_smb_gradient_anomaly_ncfile = NCDataset(vertical_smb_gradient_anomaly_filename)
    vertical_smb_gradient .= replace(vertical_smb_gradient_anomaly_ncfile["dacabfdz"][:,:,1], missing => NaN)
    
    #println("read in vertical smb gradient forcing file: " * vertical_smb_gradient_anomaly_filename)
    @info "read in vertical smb gradient forcing file: $vertical_smb_gradient_anomaly_filename"
    
    return nothing

end