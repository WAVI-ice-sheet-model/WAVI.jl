struct CoulombSlidingLaw{T <: Real, W} <: AbstractSlidingLaw 
    drag_c :: W
    reg_speed :: T
end

"""
CoulombSlidingLaw(; <kwargs>)


Keyword arguments
=================
- drag_c    : Coulomb friction coefficient
- reg_speed : regularization speed, used to prevent bed speed going to zero
"""

function CoulombSlidingLaw(; 
                        drag_c = 0.5,
                        reg_speed=1.0e-5) 
                        
    return CoulombSlidingLaw(
                            drag_c,
                            reg_speed)
end

"""
            update_β_sliding_law!(sliding_law::CoulombSlidingLaw, model::AbstractModel) 

use Coulomb sliding law to calculate basal drag
"""

function update_β_sliding_law!(sliding_law::CoulombSlidingLaw, model::AbstractModel)
    @unpack gh=model.fields
    gh.β .= (gh.drag_c .* gh.effective_pressure) ./ ( sqrt.(gh.bed_speed.^2 .+  sliding_law.reg_speed^2 ) )
    return model
end