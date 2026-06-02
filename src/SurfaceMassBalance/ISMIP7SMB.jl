export ISMIP7SMB

struct ISMIP7SMB{IC::AbstractClimateForcing, T <: Real, L <: Real} <: AbstractSurfaceMassBalance
    ISMIP7_config::IC
    reference_elevation::Union{Array{T,2}, nothing}
    vertical_smb_gradient::Union{Array{T,2}, nothing}
    smb_anomaly::Union{Array{T,2}, nothing}
    reference_smb::Union{Array{T,2}, nothing}
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
    @unpack s = model.fields.gh


    change_in_elevation = s - reference_elevation

    #set the updated smb 
    accumulation .=  reference_smb .+ smb_anomaly .+ vertical_smb_gradient .* change_in_elevation
    return nothing 


end


function update_climate_forcing!(surface_mass_balance::ISMIP7SMB, clock) 
    #do the stuff to update the forcing 
    
end