"""
InitialConditions(; 
                    initial_thickness = fill!(Array{Float64}(undef,1,1),NaN),
                    initial_grounded_fraction = fill!(Array{Float64}(undef,1,1),NaN),
                    initial_u_veloc = fill!(Array{Float64}(undef,1,1),NaN)
                    initial_v_veloc = fill!(Array{Float64}(undef,1,1),NaN)
                    initial_viscosity = fill!(Array{Float64}(undef,1,1),NaN),
                    initial_temperature = fill!(Array{Float64}(undef,1,1),NaN),
                    initial_damage = fill!(Array{Float64}(undef,1,1),NaN),
                    initial_u_veloc=fill!(Array{Float64}(undef,1,1),NaN),
                    initial_v_veloc=fill!(Array{Float64}(undef,1,1),NaN),
                    initial_basal_water_thickness = fill!(Array{Float64}(undef,1,1),NaN),
                    initial_effective_pressure = fill!(Array{Float64}(undef,1,1),NaN),
                    initial_basal_melt = fill!(Array{Float64}(undef,1,1),NaN),
                    initial_θ_ave = fill!(Array{Float64}(undef,1,1),NaN))

Construct a WAVI.jl initial conditions object. 
Unpassed arguments default to 1x1 nan matrix; unspecified initial conditions are overwritten by default values specified in Params structure when model is assembled.

Keyword arguments
=================

- 'initial_thickness': (nx x ny) matrix defining ice thickness at t = 0
- 'initial_grounded_fraction': (nx x ny) matrix defining grounded fraction at t = 0
- 'initial_u_veloc' : Initial guess for depth-averaged u velocity, including any fixed velocities (nx+1,ny). 
- 'initial_v_veloc' : Initial guess for depth-averaged v velocity, including any fixed velocities (nx,ny+1).
- 'initial_viscosity': (nx x ny x nz) matrix defining viscosity on sigma levels at t = 0
- 'initial_temperature': (nx x ny x nz) matrix defining temperature on sigma levels at t = 0
- 'initial_damage': (nx x ny x nz) matrix defining ice damage at t = 0
- 'initial_basal_water_thickness': (nx x ny) matrix defining basal water thickness at t = 0
- 'initial_effective_pressure' : (nx x ny) matrix defining effective pressure at t = 0
- 'initial_basal_melt' : (nx x ny) matrix defining basal melt rate at t = 0
- 'initial_θ_ave' : (nx x ny) matrix defining depth-averaged temperature at t = 0
"""
@with_kw struct InitialConditions{T <: Real}
    initial_thickness::Array{T,2} = fill!(Array{Float64}(undef,1,1),NaN)
    initial_grounded_fraction::Array{T,2} = fill!(Array{Float64}(undef,1,1),NaN)
    initial_u_veloc::Array{T,2} = fill!(Array{Float64}(undef,1,1),NaN)
    initial_v_veloc::Array{T,2} = fill!(Array{Float64}(undef,1,1),NaN)
    initial_viscosity::Array{T,3} = fill!(Array{Float64}(undef,1,1,1),NaN)
    initial_temperature::Array{T,3} = fill!(Array{Float64}(undef,1,1,1),NaN)
    initial_damage::Array{T,3} = fill!(Array{Float64}(undef,1,1,1),NaN)
    initial_basal_water_thickness::Array{T,2} = fill!(Array{Float64}(undef,1,1),NaN)
    initial_effective_pressure::Array{T,2} = fill!(Array{Float64}(undef,1,1),NaN)
    initial_basal_melt::Array{T,2} = fill!(Array{Float64}(undef,1,1),NaN)
    initial_θ_ave::Array{T,2} = fill!(Array{Float64}(undef,1,1),NaN)
end