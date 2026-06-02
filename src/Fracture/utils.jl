
"""
    degradation_function(d)

Degradation function used in phase field model.
"""
function degradation_function(d::T;k_reg=zero(T)) where T
  degradation = k_reg + (one(k_reg)-k_reg)*(one(d)-d)^2
  return degradation
end


positivePart(x) = 0.5*(x.+abs.(x))
negativePart(x) = 0.5*(x.-abs.(x))
degrade(x,degradation) = degradation.*positivePart(x) .+ negativePart(x)
undegrade(x,degradation) = positivePart(x)./degradation .+ negativePart(x)
