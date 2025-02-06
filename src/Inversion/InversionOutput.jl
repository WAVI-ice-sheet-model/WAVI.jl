struct InversionOutput{T <: Real}
    JKV ::  Vector{T} 
    JRMS :: Vector{T} 
end


"""
InversionOutput(; <kwargs>)

Construct a WAVI.jl parameters object for holding physical parameters.

Keyword arguments
=================
- 
"""
function InversionOutput(; JRMS = nothing,
                        JKV=nothing )
                      
  return InversionOutput(JKV,
                  JRMS
                  )
end
