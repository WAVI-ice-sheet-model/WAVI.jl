export ISMIP7Hydrofracture

using WAVI: AbstractClimateForcing
using NCDatasets

struct ISMIP7Hydrofracture{P <: String, T <: Real, CM <: Union{Array{T,2}, Nothing}} <: AbstractFracture
      hydrofracture_prefix::P
      damage_value :: T
      partially_floating_cells :: Bool
      ice_shelf_collapse_mask :: CM
      path_to_forcing :: String
end

"""
ISMIP7Hydrofracture(; <kwargs>)


Keyword arguments
=================
- hydrofracture_prefix      : prefix of the filename for the hydrofracture e.g. "ice_shelf_collapse_mask_CESM2-WACCM_ssp585_ismip7_16000m_"
- damage_value              : damage value of all floating grid cells within the ice shelf collapse mask
- partially_floating_cells  : consider partially floating cells?
- ice_shelf_collapse_mask   : mask is 0 when ice shelves are sustainable and 1 when they collapse because excessive amounts of meltwater
- path_to_forcing           : path to the forcing files folder
"""

function ISMIP7Hydrofracture(;
              hydrofracture_prefix = nothing,
              damage_value = 0.99,
              partially_floating_cells = false,
              ice_shelf_collapse_mask = nothing,
              path_to_forcing = "./"
              )

        
    #check that you've passed a hydrofracture_prefix            
    ~(hydrofracture_prefix === nothing) || throw(ArgumentError("You must pass an hydrofracture prefix"))
    #add test that the ISMIP config has the right fields?

    return ISMIP7Hydrofracture(
              hydrofracture_prefix,
              damage_value,
              partially_floating_cells,
              ice_shelf_collapse_mask, 
              path_to_forcing)
end

function update_climate_forcing!(fracture::ISMIP7Hydrofracture, grid::Grid, clock::Clock)
  @unpack dx = grid
  @unpack hydrofracture_prefix, path_to_forcing = fracture


  #get the year from clock for the forcing files
  current_time = clock.time + clock.ref_time
  current_time_string = string(Int(round(current_time)))

  # load in the ice shelf collapse mask from ISMIP7
  resolution = join([string(Int(dx)), "m"])
  ice_shelf_collapse_mask_filename = joinpath(path_to_forcing, join([hydrofracture_prefix,  current_time_string,".nc"]))
  ice_shelf_collapse_mask_ncfile   = NCDataset(ice_shelf_collapse_mask_filename)
  fracture.ice_shelf_collapse_mask .= ice_shelf_collapse_mask_ncfile["mask"][:,:,1]

  println("read in fracture mask file: " * ice_shelf_collapse_mask_filename)
  @info "read in fracture mask file: $ice_shelf_collapse_mask_filename"


  return nothing
end

function update_damage!(fracture::ISMIP7Hydrofracture,model::AbstractModel{T,N};kwargs...) where {T,N}
  @unpack gh,g3d = model.fields

  if fracture.partially_floating_cells
    g3d.Φ .= fracture.ice_shelf_collapse_mask .* (1 .- gh.grounded_fraction) .* fracture.damage_value
  else
    g3d.Φ .= fracture.ice_shelf_collapse_mask .* (gh.grounded_fraction .== 0.0) .* fracture.damage_value
  end

  return model
end

function update_strain_history!(fracture::ISMIP7Hydrofracture,model::AbstractModel{T,N};kwargs...) where {T,N}
  return model
end

function reconstruct_on_grid(fracture::ISMIP7Hydrofracture, grid::Grid)
    return ISMIP7Hydrofracture(
        fracture.hydrofracture_prefix,
        fracture.damage_value,
        fracture.partially_floating_cells,
        isnothing(fracture.ice_shelf_collapse_mask) ? zeros(grid.nx,grid.ny) :
        size(fracture.ice_shelf_collapse_mask) == (grid.nx,grid.ny) ? fracture.ice_shelf_collapse_mask :
        throw(DimensionMismatch("Size of ice shelf collapse_mask is incompatible with grid")),
        fracture.path_to_forcing)
end

function reconstruct_on_subdomain(fracture::ISMIP7Hydrofracture, grid::Grid,subdomain::NTuple{4,<: Integer})
    return ISMIP7Hydrofracture(
        fracture.hydrofracture_prefix,
        fracture.damage_value,
        fracture.partially_floating_cells,
        size(fracture.ice_shelf_collapse_mask) == size(grid)[1:2] ? fracture.ice_shelf_collapse_mask[x_start:x_end, y_start:y_end] : fracture.ice_shelf_collapse_mask,
        fracture.path_to_forcing)
end