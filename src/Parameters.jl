module Parameters

using Parameters
using WAVI.Time: compute_iterations_and_end_time

export Params, SolverParams, TimesteppingParams

# TODO: A W G? Assuming these could be functions as well as reals?
struct Params{T <: Real, A, W, G}
    dt :: T
    g :: T
    density_ice :: T
    density_ocean :: T
    gas_const :: T 
    sec_per_year  :: T
    default_thickness :: T 
    default_viscosity :: T 
    default_temperature :: T
    default_damage :: T 
    accumulation_rate :: A
    glen_a_activation_energy :: T 
    glen_a_ref :: G
    glen_temperature_ref :: T 
    glen_n :: T 
    glen_reg_strain_rate :: T 
    weertman_c :: W
    weertman_m :: T 
    weertman_reg_speed :: T 
    sea_level_wrt_geoid :: T
    minimum_thickness :: T 
    evolveShelves :: Bool
    smallHAF :: T
end


"""
    Params(; <kwargs>)

    Construct a WAVI.jl parameters object for holding physical parameters.

    Keyword arguments
    =================
    - `dt`: model timestep (NB: simulation timestep set in timestepping params, this value is updated when model embedded to the value specified in timestepping_params when passed to simulation)
    - `g`: gravitational acceleration (m^2 / s)
    - `density_ice`: ice density (kg / m^3)
    - `density_ocean`: ocean water density (kg / m^3)
    - `gas_const`: gas constant in glen b calculation
    - `sec_per_year`: seconds per year (s)
    - `default_thickness`: thickness value reverted to if no initial thickness passed (m)
    - `default_viscosity`: viscosity value reverted to if no initial thickness passed (Pa s)
    - `default_temperature`: temperature value reverted to if no initial thickness passed (K)
    - `default_damage`: damage value reverted to if no initial thickness passed (dimensionless)
    - `accumulation_rate`: uniform accumulation_rate (m/yr)
    - `glen_a_activation_energy`: activation energy in glen b calculation
    - `glen_a_ref`: array of glen a reference values used in glen b calculation
    - `glen_temperature_ref`: reference temperature using in glen b calculation
    - `glen_n`: exponent in glen b calculation
    - `glen_reg_strain_rate`: strain rate regularization value
    - `weertman_c`: basal sliding field of coefficients
    - `weertman_m`: sliding law exponent
    - `weertman_reg_speed`: regularization speed, used to prevent bed speed going to zero
    - `sea_level_wrt_geoid`: reference sea level
    - `minimum_thickness`: minimum ice thickness on model domain
    - `evolveShelves`: flag for turning on and off the evolution of the shelves in the forward run_simulation
    - `smallHAF`: small value of HAF used within update_thickness when not evolving shelves
"""
function Params(; 
    g = 9.81, 
    density_ice = 918.0,
    density_ocean = 1028.0,
    gas_const = 8.314, 
    sec_per_year =3600*24*365.25,
    default_thickness= 100.,
    default_viscosity= 1.0e7,
    default_temperature=263.15,
    default_damage= 0.0,
    accumulation_rate= 0.0,
    glen_a_activation_energy = 5.8631e+04,
    glen_a_ref= 4.9e-16 *sec_per_year * 1.0e-9,
    glen_temperature_ref= 263.15,
    glen_n = 3.0,
    glen_reg_strain_rate = 1.0e-5,
    weertman_c = 1.0e4,
    weertman_m  = 3.0,
    weertman_reg_speed = 1.0e-5,
    sea_level_wrt_geoid  = 0.0,
    minimum_thickness = 50.0,
    evolveShelves = true,
    smallHAF = 1.0)
    
    #defualt the timestep to 1.0 (will be updated when the model is embedded in a simulation)
    dt = 1.0

    return Params(
        dt, 
        g, 
        density_ice, 
        density_ocean, 
        gas_const,
        sec_per_year, 
        default_thickness, 
        default_viscosity,
        default_temperature,
        default_damage,
        accumulation_rate,
        glen_a_activation_energy,
        glen_a_ref,
        glen_temperature_ref,
        glen_n,
        glen_reg_strain_rate,
        weertman_c,
        weertman_m,
        weertman_reg_speed,
        sea_level_wrt_geoid,
        minimum_thickness,
        evolveShelves,
        smallHAF
    )
end

#structure to hold the solver parameters
@with_kw struct SolverParams{T <: Real, N <: Integer}
    n_iter_viscosity::N = 2;  @assert n_iter_viscosity ==2
    tol_picard::T = 1e-5
    maxiter_picard::N = 30
    tol_coarse::T = 1e-5
    maxiter_coarse::N = 1000
    levels::N = 3
    wavelet_threshold::T = 10.0
    nsmooth::N = 5
    smoother_omega::T = 1.0
    stencil_margin::N = 3
    super_implicitness::T = 1.0
end

struct TimesteppingParams{T <: Real, N <: Integer, TO, C, P}
            niter0 :: N      #starting iteration number
                dt :: T      #timestep
          end_time :: T      #end time of this simulation
                t0 :: T      #start time of this simulation 
        chkpt_freq :: T      #temporary checkpoint frequency
       pchkpt_freq :: T      #permanent checkpoint frequency  
       chkpt_path  :: String #path to location of permanent and tempoary checkpoint output
      n_iter_total :: TO     #total number of timesteps counting from zero
      n_iter_chkpt :: C      #number of iterations per temporary checkpoint
     n_iter_pchkpt :: P      #number of iterations per permanent checkpoint
    step_thickness :: Bool   #toggle whether to step the thickness at each timestep or not (coupling control)
end



"""
TimesteppingParams(;
                    niter0 = 0,
                    dt = 1.0,
                    end_time = 1.0,
                    n_iter_total = nothing, 
                    chkpt_freq = Inf,
                    pchkpt_freq = Inf,
                    chkpt_path = './',
                    step_thickness = true)

Construct a WAVI.jl TimesteppingParams object.
TimesteppingParams stores information relating to timestepping.

Keyword arguments
=================
- 'niter0': Iteration number of the first timestep. niter0 = 0 corresponds to a new simulation, while niter0 > 0 (positive integer) corresponds to a pickup.
- 'dt': Model timestep
- 'end_time': Simulation termination time
- 'n_iter_total': Total number of timesteps counting from zero
- 'chkpt_freq': Frequency of outputting temporary checkpoints
- 'pchkpt_freq': Frequecy with which permanent checkpoints are pass
- 'chkpt_path' : Path to location checkpoint output
- 'step_thickness': Toggle whether to update the ice thickness (true) or not (false) at each timestep
"""
function TimesteppingParams(;
                        niter0 = 0,
                        dt = 1.0,
                        end_time = nothing,
                        n_iter_total = nothing, 
                        chkpt_freq = Inf,
                        pchkpt_freq = Inf,
                        chkpt_path = "./",
                        step_thickness = true)


    #initialize t0 (really you should read start time from pickup file)
    t0 = niter0 > 0 ? niter0 * dt : 0 
    t0 = map(typeof(dt), t0)

    #check compatibility of n_iter_total and end_time, and compute them 
    end_time, n_iter_total = compute_iterations_and_end_time(end_time, n_iter_total, dt)

    #compute number of timesteps checkpoint number of timesteps
    chkpt_freq == Inf ? n_iter_chkpt = Inf : n_iter_chkpt  = round(Int, chkpt_freq/dt)
    pchkpt_freq == Inf ? n_iter_pchkpt = Inf : n_iter_pchkpt = round(Int, pchkpt_freq/dt)
    
    #check the output path ends in '/' and exists
    endswith(chkpt_path, "/") || (chkpt_path = string(chkpt_path, "/"))
    if ~isdir(chkpt_path)
        @warn string("Did not find output path ", chkpt_path, ". Any outputs will go to the working directory", pwd())
        chkpt_path = "./"
    end

    return TimesteppingParams(niter0, dt, end_time, t0, chkpt_freq, pchkpt_freq, 
                            chkpt_path,n_iter_total, n_iter_chkpt, n_iter_pchkpt, step_thickness)
end

end
