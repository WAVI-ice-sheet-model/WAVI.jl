module WAVI

#This module will export these functions and types, allowing basic use of the model.
export
    #Structures
    Model, Params, TimesteppingParams, Grid, SolverParams, InitialConditions, OutputParams, 
    Simulation, BasicParallelSpec, SharedMemorySpec,

    #Simulation controls
    update_state!, timestep!, run_simulation!,

    #Melt ratesS
    PlumeEmulator, BinfileMeltRate, UniformMeltRate, MISMIPMeltRateOne, PICO, QuadraticMeltRate, QuadraticForcedMeltRate, MeltRateExponentVariation, MeltRateExponentVariationBasins, UniformMeltUnderShelves, UniformMeltUnderShelvesBasins, 
   
    #Post-processing controls
    volume_above_floatation, height_above_floatation,

    #Abstract types
    AbstractGrid, AbstractMeltRate, AbstractParallelSpec, AbstractModel, AbstractPreconditioner


#Useful packages
using LinearAlgebra, SparseArrays, LinearMaps, Parameters,
      IterativeSolvers, Interpolations, BenchmarkTools, Reexport,
      NetCDF, JLD2, Setfield, MAT, ImageFiltering, InplaceOps,
      NonlinearSolve,SciMLNLSolve

#Import functions so they can be modified in this module.
import Base: *, size, eltype
import LinearAlgebra: ldiv!,mul!
import Setfield: @set

#Reexport Modules useful for users of the WAVI module
@reexport using JLD2
@reexport using Setfield

#Abstract types
abstract type AbstractGrid{T <: Real, N <: Integer} end
abstract type AbstractMeltRate end
abstract type AbstractParallelSpec end
abstract type AbstractModel{T <: Real, N <: Integer, M <: AbstractMeltRate, PS <: AbstractParallelSpec} end
abstract type AbstractPreconditioner{T <: Real, N <: Integer} end
#abstract type AbstractSimulation{T,N,R,A,W} end

#Type alias, just for abreviation
const MapOrMatrix{T} = Union{LinearMap{T}, AbstractMatrix{T}}

##################################################################################
#include all of the code
include("OutputParams/OutputParams.jl")
include("Grid.jl")
include("Params.jl")
include("SolverParams.jl")
include("TimesteppingParams.jl")
include("Clock.jl")
include("InitialConditions.jl")
include("KroneckerProduct.jl")
include("Wavelets/Wavelets.jl")
include("Fields/Fields.jl")
include("ParallelSpec/ParallelSpec.jl")
using .ParallelSpec
include("Models/Model.jl")
include("MeltRate/MeltRate.jl")
include("Simulations/Simulation.jl")
include("utilities.jl")


end


