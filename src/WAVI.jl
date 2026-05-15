module WAVI

# # TODO Refactor to move many of these into the submodules where they are used
using LinearAlgebra, SparseArrays, LinearMaps, Parameters,
       IterativeSolvers, Interpolations, BenchmarkTools, Reexport,
       NCDatasets, JLD2, Setfield, MAT, ImageFiltering, InplaceOps,
       NonlinearSolve, SciMLNLSolve, Enzyme

# #Import functions so they can be modified in this module.
# import Base: *, size, eltype
# import LinearAlgebra: ldiv!,mul!
# import Setfield: @set

# Abstract types
abstract type AbstractField{T <: Real, N <: Integer} end
abstract type AbstractGrid{T <: Real, N <: Integer} end
abstract type AbstractMeltRate end
abstract type AbstractFracture end
abstract type AbstractSlidingLaw end
abstract type AbstractBasalHydrology end
abstract type AbstractThermoDynamics end
abstract type AbstractSpec end
abstract type AbstractModel{T <: Real,
                            N <: Integer,
                            S <: AbstractSpec,
                            F <: AbstractField,
                            G <: AbstractGrid} end
abstract type AbstractPreconditioner{T <: Real,N <: Integer} end
abstract type AbstractSimulation end

using LinearMaps
#Type alias, just for abreviation
const MapOrMatrix{T} = Union{LinearMap{T}, AbstractMatrix{T}}

##################################################################################
#include all of the code
include("Deferred.jl")
include("Time.jl")
include("Parameters.jl")
include("KroneckerProducts.jl")
include("Grids.jl")
include("Utilities.jl")
include("Wavelets/Wavelets.jl")
include("Fields/Fields.jl")
include("SlidingLaw/SlidingLaw.jl")
include("BasalHydrology/BasalHydrology.jl")
include("ThermoDynamics/ThermoDynamics.jl")
include("MeltRates/MeltRates.jl")
include("Processes/Processes.jl")
include("Models/Models.jl")
include("Advection/Advection.jl")
include("Fracture/Fracture.jl")
include("Outputs/Outputs.jl")
include("Simulations/Simulation.jl")
include("Specs/Specs.jl")
include("Inversion/JKVsteppingParams.jl")
include("Inversion/InversionSimulation.jl")
include("Inversion/InversionOutput.jl")
include("Inversion/InversionParams.jl")
include("Inversion/DataFields.jl")
include("Inversion/Inversion.jl")

export AbstractField, AbstractGrid, AbstractMeltRate, 
  AbstractFracture, AbstractSlidingLaw , AbstractBasalHydrology,
   AbstractThermoDynamics, AbstractModel, AbstractPreconditioner,
   AbstractSpec

using .Deferred
export Collector, clear!, collect!, register_field!

using .Time
export Clock, compute_iterations_and_end_time

using .Parameters
export Params, SolverParams, TimesteppingParams

using .KroneckerProducts
export KronType, KroneckerProduct

using .Grids
export Grid

using .Utilities
export volume_above_floatation, height_above_floatation

using .Wavelets
export Preconditioner 

using .Fields
export GridField, InitialConditions

using .MeltRates
export PlumeEmulator, BinfileMeltRate, UniformMeltRate, MISMIPMeltRateOne, PICO, QuadraticMeltRate, 
QuadraticForcedMeltRate, MeltRateExponentVariation, MeltRateExponentVariationBasins, UniformMeltUnderShelves, 
UniformMeltUnderShelvesBasins

using .Processes
export update_state!, update_velocities!

using .Models
export AbstractModel, Model, update_state!

using .Outputs
export OutputParams,
    fetch_output, 
    get_spatiotemporal_var_atts, get_spatial_dimensions, get_times, get_output_as_dict, 
    make_ncfile, make_ncfile_from_filenames, 
    write_output, write_outputs, zip_output

using .Simulations
export Simulation, run_simulation!, timestep!, 
    # TODO: these probably should be in clock, processes and outputs?
    update_clock!, update_thickness!, write_vel

using .Specs
export BasicSpec, ThreadedSpec, MPISpec

# TODO These can be in their own modules as for melt rates 

#Sliding law
export WeertmanSlidingLaw, CoulombSlidingLaw, BuddSlidingLaw, TsaiSlidingLaw, TsaiBuddSlidingLaw, SchoofSlidingLaw, ZoetIversonSlidingLaw

#Basal hydrology
export NoHydrology, ConstantBasalWaterThickness, LeakyBucket, SheetOnlyGlaDS

#Thermodynamics
export NoThermoDynamics, QuadraticTemperatureApproximation, QuadraticTemperatureApproximationIcebergTest
   
#Damage models
export ConstantDamage, DruckerPragerPhaseField


end


