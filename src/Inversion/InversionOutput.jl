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
                      
    return InversionOutput(
    JKV  === nothing ? Float64[] : JKV,  # Replace nothing with an empty Float64 vector
    JRMS === nothing ? Float64[] : JRMS  # Replace nothing with an empty Float64 vector
                        )
end