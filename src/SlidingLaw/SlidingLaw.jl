#add each of the individual melt rate models
include("./coulomb.jl")

struct WeertmanSlidingLaw{T <: Real, W} <: AbstractSlidingLaw
    drag_c :: W
    weertman_m :: T
    reg_speed :: T
end

"""
WeertmanSlidingLaw(; <kwargs>)


Keyword arguments
=================
- drag_c     : Weertman friction coefficients
- weertman_m : Weertman exponent
- reg_speed  : regularization speed, used to prevent bed speed going to zero
"""

function WeertmanSlidingLaw(; 
                        drag_c = 1.0e4,
                        weertman_m  = 3.0,
                        reg_speed = 1.0e-5) 
                        
    return WeertmanSlidingLaw(
                            drag_c,
                            weertman_m,
                            reg_speed)
end

"""
            update_β_sliding_law!(sliding_law::CoulombSlidingLaw, model::AbstractModel)

use Weertman sliding law to calculate basal drag
"""

function update_β_sliding_law!(sliding_law::WeertmanSlidingLaw, model::AbstractModel)
    @unpack gh=model.fields
    gh.β .= gh.drag_c .* ( sqrt.(gh.bed_speed.^2 .+  sliding_law.reg_speed^2 ) ).^(1.0/sliding_law.weertman_m - 1.0)
    return model
end