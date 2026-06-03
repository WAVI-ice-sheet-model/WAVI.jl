export ISMIP7SMB

using WAVI: AbstractClimateForcing
using NCDatasets

struct ISMIP7SMB{IC <: AbstractClimateForcing, T <: Real, L <: Real} <: AbstractSurfaceMassBalance
    ISMIP7_config::IC
    reference_elevation::Union{Array{T,2}, Nothing}
    vertical_smb_gradient::Union{Array{T,2},Nothing}
    smb_anomaly::Union{Array{T,2}, Nothing}
    reference_smb::Union{Array{T,2}, Nothing}
end

function ISMIP7SMB(; 
                ISMIP_config = nothing,
                reference_elevation = nothing,
                vertical_smb_gradient = nothing, 
                smb_anomaly = nothing,           
                reference_smb = nothing)

    #check that you've passed an ISMIP config            
    ~(ISMIP_config === nothing) || throw(ArgumentError("You must pass an ISMIP7 config file"))
    #add test that the ISMIP config has the right fields?

    # check that you've passed reference elevation 
    ~(reference_elevation === nothing) || throw(ArgumentError("You must pass a reference elevation file"))

    #check that you've passed a vertical smb gradient
    ~(vertical_smb_gradient === nothing) || throw(ArgumentError("You must pass a reference elevation file"))

    #check that you've passed a vertical smb gradient
    ~(smb_anomaly === nothing) || throw(ArgumentError("You must pass an smb anomaly"))

    #check that you've passed a vertical smb gradient
    ~(reference_smb === nothing) || throw(ArgumentError("You must pass a reference smb"))

    return ISMIP7SMB(ISMIP_config,reference_elevation,vertical_smb_gradient, smb_anomaly,reference_smb)

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


function update_climate_forcing!(surface_mass_balance::ISMIP7SMB, grid, clock) 
    @unpack smb_anomaly = surface_mass_balance
    @unpack vertical_smb_gradient = surface_mass_balance

    #get the year from clock for the forcing files
    current_time = clock.time + clock.ref_time
    current_time_string = string(Int(round(current_time)))

    # load in the smb anomaly from ISMIP7
    resolution = join([string(Int(dx)), "m"])
    smb_anomaly_filename = join(["acabf-anomaly_AIS_", model, "_ssp", scenario, "_SDBN1-", resolution, "_v2_",  current_time_string,"_yearlyaveraged.nc"])
    smb_anomaly_ncfile   = NCDataset(smb_anomaly_filename)
    smb_anomaly .= replace(smb_anomaly_ncfile["acabf-anomaly"][:,:,:] , missing => NaN)


    # load in the vertical smb gradient from ISMIP7
    vertical_smb_gradient_anomaly_filename = join(["dacabfdz_AIS_", model, "_ssp", scenario, "_SDBN1-", resolution, "_v2_",  current_time_string,".nc"])
    vertical_smb_gradient_anomaly_ncfile = NCDataset(vertical_smb_gradient_anomaly_filename)
    vertical_smb_gradient .= replace(vertical_smb_gradient_anomaly_ncfile["dacabfdz"][:,:,:], missing => NaN)
    
    return nothing

end