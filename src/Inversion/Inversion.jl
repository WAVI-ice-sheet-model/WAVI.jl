module Inversions

export run_inversion_simulation!
export Inversion, InversionParams, JKVsteppingParams, DataFields, 
       InversionSimulation, InversionOutput

using LinearAlgebra, SparseArrays, LinearMaps, Parameters,
       IterativeSolvers, Interpolations, BenchmarkTools, 
       NCDatasets, JLD2, Setfield, MAT, ImageFiltering, InplaceOps,
       NonlinearSolve, SciMLNLSolve

using WAVI: AbstractModel, AbstractInversionSimulation
using WAVI.Outputs
import WAVI.Outputs:    zip_output
using WAVI.Time
using WAVI.Utilities
using WAVI.Processes:   update_βeff!,
                        update_βeff_on_uv_grids!,
                        update_rheological_operators!,
                        update_velocities_on_h_grid!,  
                        update_shelf_strain_rate!,
                        update_av_speed!,
                        update_bed_speed!,
                        update_surf_speed!,
                        update_surface_velocities_on_uv_grid!,
                        update_dhdt!,
                        inner_update_viscosity!,
                        update_av_viscosity!,
                        update_quadrature_falpha!,
                        update_surface_elevation!,
                        update_geometry_on_uv_grids!,
                        update_height_above_floatation!,
                        update_grounded_fraction_on_huv_grids!,
                        update_accumulation_rate!,
                        get_rhs_dirichlet!,
                        get_start_guess,
                        set_velocities!
using WAVI.Wavelets:    get_op_diag
                        

include("./InversionParams.jl")
include("./DataFields.jl")
include("./JKVsteppingParams.jl")
include("./InversionOutput.jl")
include("./InversionSimulation.jl")
include("./inversion_utilities.jl")

struct Inversion{M,IP,D,IO}
    model::M
    data::D
    inversion_params::IP
    inversion_output::IO
end

"""
        Inversion(;
            model = Model(),
            data = data,
            inversion_params = InversionParams(),
            inversion_output = InversionOutput(),
            )

Construct a WAVI.jl inversion object.

Keyword arguments
=================
    - `model`: model
    - `data`: data to be used in the inversion
    - `inversion_params`: inversion parameters
    - `inversion_output`: output from the inversion

"""
function Inversion(;
    model = Model(),
    data = DataFields(),
    inversion_params = InversionParams(),
    inversion_output = InversionOutput())

    if !all(data.ghdata.accumulation_rate .== model.params.accumulation_rate)
    println("WARNING: a different accumulation_rate is being used in model and inversion.data!")
    else
    #println("The same accumulation_rate is being used in model and inversion.data!")
    end

    inversion=Inversion(model,data,inversion_params,inversion_output)

    return inversion
end

end