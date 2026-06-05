export WeertmanSlidingLaw

struct WeertmanSlidingLaw{T <: Real, W <: Union{T,Array{T,2}}} <: AbstractSlidingLaw
    drag_coefficient :: W
    weertman_m :: T
    reg_speed :: T
end

"""
WeertmanSlidingLaw(; <kwargs>)


Keyword arguments
=================
- drag_coefficient     : Weertman friction coefficients [Pa (yr/m)^(1/weertman_m)]
- weertman_m           : Weertman exponent
- reg_speed            : regularization speed, used to prevent bed speed going to zero
"""

function WeertmanSlidingLaw(; 
                        drag_coefficient = 1.0e4,
                        weertman_m  = 3.0,
                        reg_speed = 1.0e-5) 
                        
    return WeertmanSlidingLaw(
                            drag_coefficient,
                            weertman_m,
                            reg_speed)
end

"""
            update_β_using_sliding_law!(sliding_law::WeertmanSlidingLaw, model::AbstractModel)

use Weertman sliding law to calculate basal drag coefficient
"""

function update_β_using_sliding_law!(sliding_law::WeertmanSlidingLaw, model::AbstractModel)
    @unpack gh=model.fields
    update_drag_coefficient!(model)
    gh.β .= gh.drag_coefficient .* ( sqrt.(gh.bed_speed.^2 .+  sliding_law.reg_speed^2 ) ).^(1.0/sliding_law.weertman_m - 1.0)
    return model
end

function Grids.reconstruct_on_grid(sliding_law::WeertmanSlidingLaw, grid::Grid)
    return WeertmanSlidingLaw(
        isa(sliding_law.drag_coefficient,Number) ? sliding_law.drag_coefficient*ones(grid.nx,grid.ny) : 
        size(sliding_law.drag_coefficient) == (grid.nx,grid.ny) ? sliding_law.drag_coefficient :
        throw(DimensionMismatch("Drag Coefficient does not match grid size")),
          sliding_law.weertman_m,
          sliding_law.reg_speed)
end

function Grids.reconstruct_on_subdomain(sliding_law::WeertmanSlidingLaw, grid::Grid, subdomain::NTuple{4,<: Integer})
    
    x_start,x_end,y_start,y_end = subdomain

    return WeertmanSlidingLaw(
          sliding_law.drag_coefficient[x_start:x_end, y_start:y_end],
          sliding_law.weertman_m,
          sliding_law.reg_speed)
end