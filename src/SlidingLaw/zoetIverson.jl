using WAVI.Grids: Grid


struct ZoetIversonSlidingLaw{T <: Real, C <: Union{T,Array{T,2}}, W <: Union{T,Array{T,2}}} <: AbstractSlidingLaw
    coulomb_coefficient :: C
    drag_coefficient :: W
    weertman_m :: T
    reg_speed :: T
    zoetIverson_p :: T
end

"""
ZoetIversonSlidingLaw(; <kwargs>)


Keyword arguments
=================
- coulomb_coefficient  : Coulomb friction coefficient
- drag_coefficient     : transition speed without dependence on effective pressure [m/(yr Pa)]
- weertman_m           : Weertman exponent
- reg_speed            : regularization speed, used to prevent bed speed going to zero
- zoetIverson_p        : Zoet-Iverson slip exponent
"""

function ZoetIversonSlidingLaw(;
                        coulomb_coefficient = 0.5,
                        drag_coefficient = 708.20e-6, # m/(yr Pa)
                        weertman_m  = 3.0,
                        reg_speed = 1.0e-5,
                        zoetIverson_p = 5.0) 
                        
    return ZoetIversonSlidingLaw(
                            coulomb_coefficient,
                            drag_coefficient,
                            weertman_m,
                            reg_speed,
                            zoetIverson_p)
end

"""
            update_β_using_sliding_law!(sliding_law::ZoetIversonSlidingLaw, model::AbstractModel)

use Zoet-Iverson sliding law to calculate basal drag coefficient
"""

function update_β_using_sliding_law!(sliding_law::ZoetIversonSlidingLaw, model::AbstractModel)
    @unpack gh=model.fields
    gh.β .= (sliding_law.coulomb_coefficient .* gh.effective_pressure .* ( sqrt.(gh.bed_speed.^2 .+  sliding_law.reg_speed^2 ) ).^(1.0/sliding_law.zoetIverson_p - 1.0) ./
            (( sqrt.(gh.bed_speed.^2 .+  sliding_law.reg_speed^2 ) ) .+ sliding_law.drag_coefficient .* gh.effective_pressure).^(1.0/sliding_law.zoetIverson_p))
    return model
end


function reconstruct_on_grid(sliding_law::ZoetIversonSlidingLaw, grid::Grid)
    return ZoetIversonSlidingLaw(
        isa(sliding_law.coulomb_coefficient,Number) ? sliding_law.coulomb_coefficient*ones(grid.nx,grid.ny) : sliding_law.coulomb_coefficient,
        isa(sliding_law.drag_coefficient,Number) ? sliding_law.drag_coefficient*ones(grid.nx,grid.ny) : sliding_law.drag_coefficient,
          sliding_law.weertman_m,
          sliding_law.reg_speed)
end

function reconstruct_on_subdomain(sliding_law::ZoetIversonSlidingLaw, grid::Grid, subdomain::NTuple{4,<: Integer})
    
    x_start,x_end,y_start,y_end = subdomain

    return ZoetIversonSlidingLaw(
          sliding_law.coulomb_coefficient[x_start:x_end, y_start:y_end],
          sliding_law.drag_coefficient[x_start:x_end, y_start:y_end],
          sliding_law.weertman_m,
          sliding_law.reg_speed)
end