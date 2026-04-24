struct SchoofSlidingLaw{T <: Real, C, W} <: AbstractSlidingLaw
    coulomb_coefficient :: C
    drag_coefficient :: W
    weertman_m :: T
    reg_speed :: T
end

"""
SchoofSlidingLaw(; <kwargs>)


Keyword arguments
=================
- coulomb_coefficient  : Cmax = Iken's bound = tan(maximum up-slope angle of the bed in flow direction)
- drag_coefficient     : Schoof friction coefficients [Pa (yr/m)^(1/3)]
- weertman_m           : Weertman exponent
- reg_speed            : regularization speed, used to prevent bed speed going to zero
"""

function SchoofSlidingLaw(;
                        coulomb_coefficient = 0.2,
                        drag_coefficient = 3.165e6, # 3.165e6 = Pa (yr/m)^(1/3) = 1.0e3 MPa (s/m)^(1/3)
                        weertman_m  = 3.0,
                        reg_speed = 1.0e-5) 
                        
    return SchoofSlidingLaw(
                            coulomb_coefficient,
                            drag_coefficient,
                            weertman_m,
                            reg_speed)
end

"""
            update_β_using_sliding_law!(sliding_law::SchoofSlidingLaw, model::AbstractModel)

use Schoof sliding law to calculate basal drag coefficient
"""

function update_β_using_sliding_law!(sliding_law::SchoofSlidingLaw, model::AbstractModel)
    @unpack gh=model.fields
    gh.β .= (sliding_law.drag_coefficient .* ( sqrt.(gh.bed_speed.^2 .+  sliding_law.reg_speed^2 ) ).^(1.0/sliding_law.weertman_m - 1.0) ./
            (1 .+ (sliding_law.drag_coefficient ./ (sliding_law.coulomb_coefficient .* gh.effective_pressure)).^(sliding_law.weertman_m) .* ( sqrt.(gh.bed_speed.^2 .+  sliding_law.reg_speed^2 ) )).^(1.0/sliding_law.weertman_m))
    return model
end