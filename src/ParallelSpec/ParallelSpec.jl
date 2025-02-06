module ParallelSpec

export BasicParallelSpec

using Parameters, Setfield

using WAVI: AbstractParallelSpec, AbstractModel

struct BasicParallelSpec <: AbstractParallelSpec end

include("shared_memory_spec.jl")

end