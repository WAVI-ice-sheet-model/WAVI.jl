export ISMIP7hydrofracture

struct ISMIP7hydrofracture{T <: Real, CM <: Union{Array{T,2}, Nothing}} <: AbstractFracture
      damage_value :: T
      partially_floating_cells :: Bool
      ice_shelf_collapse_mask :: CM
      ice_shelf_collapse_file :: String
end

"""
ISMIP7hydrofracture(; <kwargs>)


Keyword arguments
=================
- damage_value              : damage value of all floating grid cells within the ice shelf collapse mask
- partially_floating_cells  : consider partially floating cells?
- ice_shelf_collapse_mask   : mask is 0 when ice shelves are sustainable and 1 when they collapse because excessive amounts of meltwater
- ice_shelf_collapse_file   : NetCDF file containing the ice shelf collapse mask
"""

function ISMIP7hydrofracture(;
              damage_value = 0.99,
              partially_floating_cells = false,
              ice_shelf_collapse_mask = nothing,
              ice_shelf_collapse_file = nothing
              )

    return ISMIP7hydrofracture(
              damage_value,
              partially_floating_cells,
              ice_shelf_collapse_mask,
              ice_shelf_collapse_file)
end

function update_climate_forcing!(fracture::ISMIP7hydrofracture, clock)
  # ice_shelf_collapse_file = "/data/icesheet_output/miradh/ISMIP7/ISMIP7/AIS/CESM2-WACCM/ssp585/fracture/ice_shelf_collapse_mask_cesm2waccm_ssp585_ismip7_8km.nc"
  ds = NCDataset(fracture.ice_shelf_collapse_file)
  timeMask = ds["time"][:]
  timeInd = findfirst(timeMask .== Int(round(clock.ref_time + clock.time)))
  fracture.ice_shelf_collapse_mask .= ds["mask"][:,:,timeInd]
  close(ds)

  return model
end

function update_damage!(fracture::ISMIP7hydrofracture,model::AbstractModel{T,N};kwargs...) where {T,N}
  @unpack gh,g3d = model.fields

  if partially_floating_cells
    g3d.Φ .= fracture.ice_shelf_collapse_mask .* (1 .- gh.grounded_fraction) .* fracture.damage_value
  else
    g3d.Φ .= fracture.ice_shelf_collapse_mask .* (gh.grounded_fraction .== 0.0) .* fracture.damage_value
  end

  return model
end

function update_strain_history!(fracture::ISMIP7hydrofracture,model::AbstractModel{T,N};kwargs...) where {T,N}
  return model
end