export CoulombSlidingLaw
using WAVI.Grids: Grid


struct CoulombSlidingLaw{T <: Real, C <: Union{T,Array{T,2}}} <: AbstractSlidingLaw 
    coulomb_coefficient :: C
    reg_speed :: T
end

"""
CoulombSlidingLaw(; <kwargs>)


Keyword arguments
=================
- coulomb_coefficient    : Coulomb friction coefficient
- reg_speed : regularization speed, used to prevent bed speed going to zero
"""

function CoulombSlidingLaw(; 
                        coulomb_coefficient = 0.5,
                        reg_speed=1.0e-5) 
                        
    return CoulombSlidingLaw(
                            coulomb_coefficient,
                            reg_speed)
end

"""
            update_β_using_sliding_law!(sliding_law::CoulombSlidingLaw, model::AbstractModel) 

use Coulomb sliding law to calculate basal drag coefficient
"""

function update_β_using_sliding_law!(sliding_law::CoulombSlidingLaw, model::AbstractModel)
    @unpack gh=model.fields
    gh.β .= (sliding_law.coulomb_coefficient .* gh.effective_pressure) ./ ( sqrt.(gh.bed_speed.^2 .+  sliding_law.reg_speed^2 ) )
    return model
end


function reconstruct_on_grid(sliding_law::CoulombSlidingLaw, grid::Grid)
    return CoulombSlidingLaw(
        isa(sliding_law.coulomb_coefficient,Number) ? sliding_law.coulomb_coefficient*ones(grid.nx,grid.ny) : sliding_law.coulomb_coefficient,
          sliding_law.reg_speed)
end

function reconstruct_on_subdomain(sliding_law::CoulombSlidingLaw, grid::Grid, subdomain::NTuple{4,<: Integer})
    
    x_start,x_end,y_start,y_end = subdomain

    return CoulombSlidingLaw(
          sliding_law.coulomb_coefficient[x_start:x_end, y_start:y_end],
          sliding_law.reg_speed)
end