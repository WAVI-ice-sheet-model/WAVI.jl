module Specs

using WAVI: AbstractSpec

export AbstractDecompSpec

abstract type AbstractDecompSpec <: AbstractSpec end

include("basic.jl")
#include("threaded.jl")
#include("mpi.jl")

end