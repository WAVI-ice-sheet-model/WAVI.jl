export TsaiBuddSlidingLaw

struct TsaiBuddSlidingLaw{T <: Real, C <: Union{T,Array{T,2}}, W <: Union{T,Array{T,2}}} <: AbstractSlidingLaw
    coulomb_coefficient :: C
    drag_coefficient :: W
    weertman_m :: T
    reg_speed :: T
    budd_q :: T
end

"""
TsaiBuddSlidingLaw(; <kwargs>)


Keyword arguments
=================
- coulomb_coefficient  : Coulomb friction coefficient
- drag_coefficient     : Budd friction coefficients [(yr/m)^(1/weertman_m)]
- weertman_m           : Weertman exponent
- reg_speed            : regularization speed, used to prevent bed speed going to zero
- budd_q               : Budd exponent
"""

function TsaiBuddSlidingLaw(; 
                        coulomb_coefficient = 0.5,
                        drag_coefficient = 0.117, # (yr/m)^(1/3) = 37.01 (s/m)^(1/3)
                        weertman_m  = 3.0,
                        reg_speed = 1.0e-5,
                        budd_q = 1.0) 
                        
    return TsaiBuddSlidingLaw(
                            coulomb_coefficient,
                            drag_coefficient,
                            weertman_m,
                            reg_speed,
                            budd_q)
end

"""
            update_β_using_sliding_law!(sliding_law::TsaiBuddSlidingLaw, model::AbstractModel)

use Tsai-Budd sliding law to calculate basal drag coefficient
"""

function update_β_using_sliding_law!(sliding_law::TsaiBuddSlidingLaw, model::AbstractModel)
    @unpack gh=model.fields
    gh.β .= (min.((sliding_law.coulomb_coefficient .* gh.effective_pressure) ./ ( sqrt.(gh.bed_speed.^2 .+  sliding_law.reg_speed^2 ) ), 
                sliding_law.drag_coefficient .* ( sqrt.(gh.bed_speed.^2 .+  sliding_law.reg_speed^2 ) ).^(1.0/sliding_law.weertman_m - 1.0) .* gh.effective_pressure.^(sliding_law.budd_q)))
    return model
end


function reconstruct_on_grid(sliding_law::TsaiBuddSlidingLaw, grid::Grid)
    return TsaiBuddSlidingLaw(
        isa(sliding_law.coulomb_coefficient,Number) ? sliding_law.coulomb_coefficient*ones(grid.nx,grid.ny) : 
        size(sliding_law.coulomb_coefficient) == (grid.nx,grid.ny) ? sliding_law.coulomb_coefficient :
        throw(DimensionMismatch("Coulomb Coefficient does not match grid size")),
        isa(sliding_law.drag_coefficient,Number) ? sliding_law.drag_coefficient*ones(grid.nx,grid.ny) : 
        size(sliding_law.drag_coefficient) == (grid.nx,grid.ny) ? sliding_law.drag_coefficient :
        throw(DimensionMismatch("Drag Coefficient does not match grid size")),
        sliding_law.weertman_m,
        sliding_law.reg_speed,
        budd_q)
end

function reconstruct_on_subdomain(sliding_law::TsaiBuddSlidingLaw, grid::Grid, subdomain::NTuple{4,<: Integer})
    
    x_start,x_end,y_start,y_end = subdomain

    return TsaiBuddSlidingLaw(
          sliding_law.coulomb_coefficient[x_start:x_end, y_start:y_end],
          sliding_law.drag_coefficient[x_start:x_end, y_start:y_end],
          sliding_law.weertman_m,
          sliding_law.reg_speed,
          budd_q)
end