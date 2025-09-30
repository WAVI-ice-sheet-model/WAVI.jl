struct InversionOutput{T <: Real}
    JKV ::  Vector{T} 
    JRMS :: Vector{T} 
end


"""
InversionOutput(; <kwargs>)

Construct a WAVI.jl InversionOutput object for holding the inversion output parameters.


Keyword arguments
=================
- JKV:                   Vector of JKV values
- JRMS:                  Vector of JRMS values
"""
function InversionOutput(; JRMS = nothing,
                        JKV=nothing )
                      
 return InversionOutput(
    JKV  === nothing ? Float64[] : JKV,  # Replace nothing with an empty Float64 vector
    JRMS === nothing ? Float64[] : JRMS  # Replace nothing with an empty Float64 vector
                        )
end