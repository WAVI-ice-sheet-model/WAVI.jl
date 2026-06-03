export ISMIP7MeltRate

using WAVI: AbstractClimateForcing
using NCDatasets 

struct ISMIP7MeltRate{IC <: AbstractClimateForcing, T <: Real} <: AbstractMeltRate
    ISMIP7_config::IC                   #ISMIP 7 config file, specifies the GCM and scenario
    K::T                                #tuning parameter (see here for info: https://drive.google.com/file/d/1SygMQte-7XgKj4e-Hpyj1tHi6XguShCP/view and https://tc.copernicus.org/articles/16/4931/2022/tc-16-4931-2022.pdf)
    shelf_slope::T                      #average slope of ice shelf base (can be local in ISMIP7, but we're doing with the average)
    ρ_ocean::T                          #density of the ocean
    ρ_ice::T                            #density of the ice 
    c_ocean::T                             #thermal heat paracity of the ocean 
    L_ice::T                            #ice latent heat of fusion 
    β_s::T                              #Salt contraction coefficient
    g::T                                #gravity 
    f::T                                #Coriolis parameter
    S_loc::Union{Array{T,3}, Nothing}   #array to hold the 3d salinity, to be read in from ISMIP7 forcing
    T_loc::Union{Array{T,3}, Nothing}   #array to hold the 3d temperature, to be read in from ISMIP7 forcing
    Tf_loc::Union{Array{T,3}, Nothing}  #array to hold the 3d thermal forcing, to be read in from ISMIP7 forcing
    z_forcing::Union{Array{T,1}, Nothing} #array to hold the z co-ordinates from the ISMIP7 forcing
    melt_partial_cell::Bool             #Flag for melting applied to partial cells or not
end

"""
    ISMIP7MeltRate(;, kwargs)

Construct an ISMIP7MeltRate object to prescribe the melt rate in WAVI

Keyword arguments
=================
- ISMIP7_config: specifies the ISMIP 7 config (GCM and scenario)
- K: melt rate tuning parameter (see here for info: https://drive.google.com/file/d/1SygMQte-7XgKj4e-Hpyj1tHi6XguShCP/view and https://tc.copernicus.org/articles/16/4931/2022/tc-16-4931-2022.pdf)
- shelf_slope: average slope of ice shelf base (can be local in ISMIP7, but we're doing with the average)
- ρ_ocean: density of the ocean water
- ρ_ice: density of the ice 
- c_ocean: thermal heat paracity of the ocean 
- L_ice: ice latent heat of fusion 
- Salt: contraction coefficient
- g: gravity 
- f: Coriolis parameter
- S_loc: array to hold the 3d salinity, to be read in from ISMIP7 forcing
- array: to hold the 3d temperature, to be read in from ISMIP7 forcing
- Tf_loc: array to hold the 3d thermal forcing, to be read in from ISMIP7 forcing
- z_forcing: array to hold the z co-ordinates from the ISMIP7 forcing
- melt_partial_cell: Flag for melting applied to partial cells or not
"""
function ISMIP7MeltRate(; 
                        ISMIP7_config = nothing,
                        K = 1.0e-5,
                        shelf_slope = 1.e-3,
                        ρ_ocean = 1028.0,
                        ρ_ice = 918.0,
                        c_ocean = 3974.0,
                        L_ice = 3.34e5,
                        β_s = 7.7e-4,
                        g = 9.81,
                        f = -1.0e-4,
                        S_loc = nothing,
                        T_loc = nothing,
                        Tf_loc = nothing, 
                        z_forcing = nothing, 
                        melt_partial_cell = false)

    #check that you've passed an ISMIP config            
    ~(ISMIP_config === nothing) || throw(ArgumentError("You must pass an ISMIP7 config file"))

    return   ISMIP7MeltRate(ISMIP7_config, K, shelf_slope,ρ_ocean, ρ_ice,c_ocean, L_ice,β_s, g,f,S_loc, T_loc, Tf_loc, z_forcing, melt_partial_cell)
end


function update_shelf_melt_rate!(ISMIP7_melt_rate::ISMIP7MeltRate, fields, grid, clock)

    @unpack b, shelf_basal_melt, grounded_fraction, h = fields.gh 
    @unpack K, shelf_slope, ρ_ocean, ρ_ice, c_ocean, L_ice, β_s, g, f, S_loc, T_loc, Tf_loc, z_forcing, melt_partial_cell = ISMIP7_melt_rate
    
    #compute the ice draft
    zb = b .* (grounded_fraction .== 1) + - ρ_ice / ρ_w .* h .* (grounded_fraction .< 1)

    #compute the local thermal forcing and salinity from the climate forcing and ice draft
    idx = [argmin(abs.(z_forcing .- d)) for d in b] #gives an nx * ny array indices which are closest depth to z in the forcing.
    Tf_local_shelf = [Tf_loc[i,j,idx[i,j]] for i in axes(Tf_loc,1), j in axes(Tf_loc,2)] #fills the array of thermal forcing 
    S_local_shelf = [S_loc[i,j,idx[i,j]] for i in axes(S_loc,1), j in axes(S_loc,2)] #fills the array of thermal forcing 


    #set the shelf melt rate
    secs_per_year = 365.25*24*60^2
    if melt_partial_cell
        shelf_basal_melt[:] .= K * secs_per_year * shelf_slope * ρ_ocean / ρ_ice * (c_ocean/L_ice)^2 * β_s * 2 * g /abs(f) * S_local_shelf[:] .* abs.(Tf_local_shelf[:]) .* Tf_local_shelf[:] .*  (1 .- grounded_fraction)

    elseif ~(melt_partial_cell)
        shelf_basal_melt[grounded_fraction .== 0] .= K * secs_per_year * shelf_slope * ρ_ocean / ρ_ice * (c_ocean/L_ice)^2 * β_s * 2 * g /abs(f) * S_local_shelf[grounded_fraction .== 0] .* abs.(Tf_local_shelf[grounded_fraction .== 0]) .* Tf_local_shelf[grounded_fraction .== 0]
        shelf_basal_melt[.~(grounded_fraction .== 0)] .= 0
    end


    return nothing
end

function update_climate_forcing!(ISMIP7_melt_rate::ISMIP7MeltRate, grid, clock)
    @unpack S_loc, T_loc, Tf_loc, z_forcing = ISMIP7_melt_rate
    @unpack dx = grid

    #get the year from clock for the forcing files
    current_time = clock.time + clock.ref_time
    current_time_string = string(Int(round(current_time)))

    # load in the salinity from ISMIP7
    resolution = join([string(Int(dx)), "m"])
    salinity_filename = join(["so_AIS_", model, "_ssp", scenario, "ocean_", resolution, "_v3_",  current_time_string,"_yearlyaveraged.nc"])
    salinity_ncfile   = NCDataset(salinity_filename)
    S_loc .= replace(salinity_ncfile["so"][:,:,:] , missing => NaN)

    # load in the temperature from ISMIP7
    temperature_filename = join(["thetao_AIS_", model, "_ssp", scenario, "ocean_", resolution, "_v3_",  current_time_string,"_yearlyaveraged.nc"])
    temperature_ncfile   = NCDataset(temperature_filename)
    T_loc .= replace(temperature_ncfile["thetao"][:,:,:] , missing => NaN)

    # load in the temperature from ISMIP7
    thermal_forcing_filename = join(["so_AIS_", model, "_ssp", scenario, "ocean_", resolution, "_v3_",  current_time_string,"_yearlyaveraged.nc"])
    thermal_forcing_ncfile   = NCDataset(thermal_forcing_filename)
    Tf_loc .= replace(thermal_forcing_ncfile["tf"][:,:,:] , missing => NaN)

    # read in the z co-ordinates 
    z_forcing .= replace(thermal_forcing_ncfile["z"][:] , missing => NaN)

end