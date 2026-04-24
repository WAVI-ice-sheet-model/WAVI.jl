struct TsaiSlidingLaw{T <: Real, C, W} <: AbstractSlidingLaw
    coulomb_coefficient :: C
    drag_coefficient :: W
    weertman_m :: T
    reg_speed :: T
end

"""
TsaiSlidingLaw(; <kwargs>)


Keyword arguments
=================
- coulomb_coefficient  : Coulomb friction coefficient
- drag_coefficient     : Weertman friction coefficients [Pa (yr/m)^(1/weertman_m)]
- weertman_m           : Weertman exponent
- reg_speed            : regularization speed, used to prevent bed speed going to zero
"""

function TsaiSlidingLaw(; 
                        coulomb_coefficient = 0.5,
                        drag_coefficient = 1.0e4,
                        weertman_m  = 3.0,
                        reg_speed = 1.0e-5) 
                        
    return TsaiSlidingLaw(
                            coulomb_coefficient,
                            drag_coefficient,
                            weertman_m,
                            reg_speed)
end

"""
            update_β_using_sliding_law!(sliding_law::TsaiSlidingLaw, model::AbstractModel)

use Tsai sliding law to calculate basal drag coefficient
"""

function update_β_using_sliding_law!(sliding_law::TsaiSlidingLaw, model::AbstractModel)
    @unpack gh=model.fields
    update_drag_coefficient!(model)
    gh.β .= (min.((sliding_law.coulomb_coefficient .* gh.effective_pressure) ./ ( sqrt.(gh.bed_speed.^2 .+  sliding_law.reg_speed^2 ) ), 
                gh.drag_coefficient .* ( sqrt.(gh.bed_speed.^2 .+  sliding_law.reg_speed^2 ) ).^(1.0/sliding_law.weertman_m - 1.0)))
    return model
end