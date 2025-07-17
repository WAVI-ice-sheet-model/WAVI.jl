struct UniformMeltUnderShelves{T <: Real} <: AbstractMeltRate 
    melt_constant  :: T
    melt_partial_cell :: Bool
    ρi :: T
    ρw :: T
end

"""
    function UniformMeltUnderShelves(; <kwargs>)


Keyword arguments
=================
- melt_constant : uniform melt rate (defaults to zero)
- melt_partial_cell : whether to melt under partially grounded cells or not (default false)
- ρi : Ice density
- ρw : Water density
"""
UniformMeltUnderShelves(; melt_constant = 0.0, melt_partial_cell= false, ρi = 918.0, ρw = 1028.0) = UniformMeltUnderShelves(melt_constant,melt_partial_cell,ρi, ρw)


"""
    update_shelf_melt_rate(shelf_melt_rate::UniformMeltUnderShelves, fields, grid, clock) 

Update the melt rate under ice shelves for the UniformMeltUnderShelves type
"""
function update_shelf_melt_rate!(shelf_melt_rate::UniformMeltUnderShelves, fields, grid, clock) 
    m = zeros(grid.nx,grid.ny);

    if (shelf_melt_rate.melt_partial_cell)  #partial cell melting 
        m[:] .=   shelf_melt_rate.melt_constant.* (1 .- fields.gh.grounded_fraction[:])
        
    elseif ~shelf_melt_rate.melt_partial_cell #no partial cell melting
        m[:] .=  shelf_melt_rate.melt_constant.* (1 .- fields.gh.grounded_fraction[:])
        m[.~(fields.gh.grounded_fraction .== 0)] .= 0 
    end
    
    fields.gh.shelf_basal_melt[:] .= m[:]
end


