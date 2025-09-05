export UniformMeltUnderShelvesBasins

struct UniformMeltUnderShelvesBasins{T <: Real} <: AbstractMeltRate 
    melt_constant_basin_1  :: T
    melt_constant_basin_2  :: T
    basinID_1 :: T                 # basin ID value for basin 1
    basinID_2 :: T                 # basin ID value for basin 2
    melt_partial_cell :: Bool
    ρi :: T
    ρw :: T
end

"""
    function UniformMeltUnderShelvesBasins(; <kwargs>)

Construct a melt rate to specify the melt rate as a constant under floating ice and zero elsewhere

Keyword arguments
=================
- melt_constant_basin_1 : uniform melt rate in basin 1 (defaults to zero)
- melt_constant_basin_2 : uniform melt rate in basin 2  (defaults to zero)
- basinID_1: ID number of basin 1 
- basinID_2: ID number of basin 2 
- melt_partial_cell : whether to melt under partially grounded cells or not (default false)
- ρi : Ice density
- ρw : Water density
"""
UniformMeltUnderShelvesBasins(; melt_constant_basin_1 = 0.0, melt_constant_basin_2 = 0.0,  basinID_1=1.0, basinID_2=2.0, melt_partial_cell= false, ρi = 918.0, ρw = 1028.0) = UniformMeltUnderShelvesBasins(melt_constant_basin_1,melt_constant_basin_2, basinID_1, basinID_2, melt_partial_cell,ρi, ρw)

"""
    update_melt_rate(melt_rate::UniformMeltUnderShelvesBasins, fields, grid) 

Update the melt rate when for the UniformMeltUnderShelvesBasins type
"""
function update_melt_rate!(melt_rate::UniformMeltUnderShelvesBasins, fields, grid, clock) 

  m = zeros(grid.nx,grid.ny);

    if (melt_rate.melt_partial_cell)  #partial cell melting 
        m[(grid.basin_ID .==melt_rate.basinID_1)] .=   melt_rate.melt_constant_basin_1.* (1 .- fields.gh.grounded_fraction[(grid.basin_ID .==melt_rate.basinID_1)])
        m[(grid.basin_ID .==melt_rate.basinID_2)] .=   melt_rate.melt_constant_basin_2.* (1 .- fields.gh.grounded_fraction[(grid.basin_ID .==melt_rate.basinID_2)])
    elseif ~melt_rate.melt_partial_cell #no partial cell melting
        m[(grid.basin_ID .==melt_rate.basinID_1)] .=  melt_rate.melt_constant_basin_1.* (1 .- fields.gh.grounded_fraction[(grid.basin_ID .==melt_rate.basinID_1)])
        m[(grid.basin_ID .==melt_rate.basinID_2)] .=  melt_rate.melt_constant_basin_2.* (1 .- fields.gh.grounded_fraction[(grid.basin_ID .==melt_rate.basinID_2)])
        m[.~(fields.gh.grounded_fraction .== 0)] .= 0 
    end

    fields.gh.basal_melt[:] .= m[:]
end

