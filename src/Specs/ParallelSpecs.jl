module ParallelSpecs

using MiniWAVI: AbstractSpec

export AbstractDecompSpec

abstract type AbstractDecompSpec <: AbstractSpec end

include("basic.jl")
include("mpi.jl")

end