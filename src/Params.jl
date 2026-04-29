struct Params{T <: Real, A, G, H, E}
                      dt :: T
                       g :: T
             density_ice :: T
      density_freshwater :: T
           density_ocean :: T
               gas_const :: T 
           sec_per_year  :: T
       default_thickness :: T 
       default_viscosity :: T 
     default_temperature :: T
          default_damage :: T 
  default_strain_history :: T
    ice_tensile_strength :: T
ice_compressive_strength :: T
 critical_elastic_energy :: T
      phase_field_length :: T
     energy_release_rate :: T
degradation_regularisation :: T
       accumulation_rate :: A
glen_a_activation_energy :: T 
              glen_a_ref :: G
    glen_temperature_ref :: T 
                  glen_n :: T 
    glen_reg_strain_rate :: T 
          elastic_lambda :: T 
              elastic_mu :: T 
     sea_level_wrt_geoid :: T
       minimum_thickness :: T 
           evolveShelves :: Bool
                smallHAF :: T
   basal_water_thickness :: T
   hydraulic_potential_b :: H
      effective_pressure :: E
              basal_melt :: T
 default_temperature_ave :: T
end


"""
Params(; <kwargs>)

Construct a WAVI.jl parameters object for holding physical parameters.

Keyword arguments
=================
- `dt`: model timestep (NB: simulation timestep set in timestepping params, this value is updated when model embedded to the value specified in timestepping_params when passed to simulation)
- `g`: gravitational acceleration (m / s^2)
- `density_ice`: ice density (kg / m^3)
- `density_freshwater`: freshwater density (kg / m^3)
- `density_ocean`: ocean water density (kg / m^3)
- `gas_const`: gas constant in glen b calculation
- `sec_per_year`: seconds per year (s)
- `default_thickness`: thickness value reverted to if no initial thickness passed (m)
- `default_viscosity`: viscosity value reverted to if no initial thickness passed (Pa s)
- `default_temperature`: temperature value reverted to if no initial thickness passed (K)
- `default_damage`: damage value reverted to if no initial damage passed (dimensionless)
- `default_strain_history`: maximum previous strain energy reverted to if no initial value passed (Pa)
- `ice_tensile_strength`: uniaxial tensile strength of ice (Pa)
- `ice_compressive_strength`: uniaxial conpressive strength of ice (Pa)
- `critical_elastic_energy`: critical elastic strain energy for fracture nucleation (Pa)
- `phase_field_length`: phase field length for fracture (m)
- `energy_release_rate`: energy release rate for fracture (N/m)
- `degradation_regularisation`: smallest degradation allowed
- `accumulation_rate`: uniform accumulation_rate (m/yr)
- `glen_a_activation_energy`: activation energy in glen b calculation
- `glen_a_ref`: array of glen a reference values used in glen b calculation
- `glen_temperature_ref`: reference temperature using in glen b calculation
- `glen_n`: exponent in glen b calculation
- `glen_reg_strain_rate`: strain rate regularization value
- `elastic_lambda`: First elastic Lame parameter (Pa)
- `elastic_mu`: Second elastic Lame parameter (Pa) 
- `sea_level_wrt_geoid`: reference sea level
- `minimum_thickness`: minimum ice thickness on model domain
- `evolveShelves`: flag for turning on and off the evolution of the shelves in the forward run_simulation
- `smallHAF`: small value of HAF used within update_thickness when not evolving shelves
- `basal_water_thickness` : basal water thickness (m)
- `hydraulic_potential_b` : hydraulic_potential at the bed (Pa)
- `effective_pressure` : effective pressure (Pa)
- `basal_melt` : basal melt rate (m/yr)
- `default_temperature_ave` : depth-averaged temperature (K)
"""
function Params(; g = 9.81, 
                  density_ice = 918.0,
                  density_freshwater = 1000.0,
                  density_ocean = 1028.0,
                  gas_const = 8.314, 
                  sec_per_year =3600*24*365.25,
                  default_thickness= 100.,
                  default_viscosity= 1.0e7,
                  default_temperature=263.15,
                  default_damage= 0.0,
                  default_strain_history = 0.0,
                  ice_tensile_strength = 262.9e3,
                  ice_compressive_strength =  496.7e3,
                  critical_elastic_energy = 1.0,
                  phase_field_length = 1.0,
                  energy_release_rate = 1.0,
                  degradation_regularisation = 1e-2,
                  accumulation_rate= 0.0,
                  glen_a_activation_energy = 5.8631e+04,
                  glen_a_ref= 4.9e-16 *sec_per_year * 1.0e-9,
                  glen_temperature_ref= 263.15,
                  glen_n = 3.0,
                  glen_reg_strain_rate = 1.0e-5,
                  elastic_lambda = 6.54e9,
                  elastic_mu = 3.52e9,
                  sea_level_wrt_geoid  = 0.0,
                  minimum_thickness = 50.0,
                  evolveShelves = true,
                  smallHAF = 1.0,
                  basal_water_thickness = 0.0,
                  hydraulic_potential_b = 0.0,
                  effective_pressure = 1.0e6,
                  basal_melt = 0.0,
                  default_temperature_ave = 253.15)
                      
  #default the timestep to 1.0 (will be updated when the model is embedded in a simulation)
  dt = 1.0

  return Params(
                  dt, 
                  g, 
                  density_ice,
                  density_freshwater, 
                  density_ocean, 
                  gas_const,
                  sec_per_year, 
                  default_thickness, 
                  default_viscosity,
                  default_temperature,
                  default_damage,
                  default_strain_history,
                  ice_tensile_strength,
                  ice_compressive_strength,
                  critical_elastic_energy,
                  phase_field_length,
                  energy_release_rate,
                  degradation_regularisation,
                  accumulation_rate,
                  glen_a_activation_energy,
                  glen_a_ref,
                  glen_temperature_ref,
                  glen_n,
                  glen_reg_strain_rate,
                  elastic_lambda,
                  elastic_mu,
                  sea_level_wrt_geoid,
                  minimum_thickness,
                  evolveShelves,
                  smallHAF,
                  basal_water_thickness,
                  hydraulic_potential_b,
                  effective_pressure,
                  basal_melt,
                  default_temperature_ave
                  )
end
