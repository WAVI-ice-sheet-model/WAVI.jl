module WAVI

#This module will export these functions and types, allowing basic use of the model.
export
    #Structures
    TimesteppingParams, OutputParams, Simulation,

    #Simulation controls
    timestep!, run_simulation!,

    #Melt ratesS
    PlumeEmulator, BinfileMeltRate, UniformMeltRate, MISMIPMeltRateOne, PICO, QuadraticMeltRate, QuadraticForcedMeltRate, MeltRateExponentVariation, MeltRateExponentVariationBasins, UniformMeltUnderShelves, UniformMeltUnderShelvesBasins, 
   
    #Post-processing controls
    volume_above_floatation, height_above_floatation,

    #Abstract types
    AbstractField, AbstractGrid, AbstractMeltRate, AbstractModel, AbstractPreconditioner, AbstractSpec


# Abstract types
abstract type AbstractField{T <: Real} end
abstract type AbstractGrid{T <: Real} end
abstract type AbstractMeltRate end
abstract type AbstractModel{T <: Real, 
                            N <: Integer, 
                            S <: AbstractSpec,
                            F <: AbstractField,
                            G <: AbstractGrid,
                            M <: AbstractMeltRate} end
abstract type AbstractPreconditioner{T <: Real} end
abstract type AbstractSpec end

#Type alias, just for abreviation
const MapOrMatrix{T} = Union{LinearMap{T}, AbstractMatrix{T}}

##################################################################################
#include all of the code
include("Time.jl")
include("Parameters.jl")
include("KroneckerProducts.jl")
include("Grids.jl")
include("Utilities.jl")
include("Fields/Fields.jl")
include("Processes/Processes.jl")
include("Models/Models.jl")
include("Specs/Specs.jl")
include("OutputParams/OutputParams.jl")
include("TimesteppingParams.jl")
include("Wavelets/Wavelets.jl")
include("MeltRates/MeltRates.jl")
include("Simulations/Simulation.jl")

using .Parameters
export Params, SolverParams

using .Grids
export Grid

using .Fields
export InitialConditions

using .Processes
export update_state!, update_velocities!

using .Models
export AbstractModel, Model,
    update_state!

using .ParallelSpecs
export BasicSpec, ThreadedSpec, MPISpec

end


