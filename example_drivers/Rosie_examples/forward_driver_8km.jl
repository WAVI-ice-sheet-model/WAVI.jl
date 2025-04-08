#This is a test 8km relaxation script: 

using WAVI 

function forward_driver_8km()

#Load up the model from the inversion:
n_ref_relxation=10;
  
file_path="/data/glacio/icesheets/shelf_thin/wip/chll1/WAVI_inversion/WAVI.jl/outputs_8km_relax_test/"
#PChkpt_0000000015.jld2"

    n_ref_string =  lpad(n_ref_relxation, 10, "0"); #filename as a string with 10 digits
        println("Requested relaxation is iteration $n_ref_relxation. Looking for relaxation file...")
    #    try 
            #for running on vscode vis workstation:
            filename = string(file_path, "/PChkpt_", n_ref_string, ".jld2")
            #for running with ensembler:
            #filename = string("PChkpt_", n_iter_string, ".jld2")
            println("Looking for file: ", filename)

            if isfile(filename)
                println("relaxation file exists")
                sim_load = load(filename, "simulation")
                println("relaxation pickup successful")
            else
                error("File $filename not found!")
            end
            

#load up grid used for inversion:  
grid = sim_load.model.grid
          
# set up initial conditions from inversion:        
initial_conditions = InitialConditions(initial_thickness =  sim_load.model.fields.gh.h,
                                        initial_viscosity = sim_load.model.fields.g3d.η,
                                        initial_temperature = sim_load.model.fields.g3d.θ,
                                        initial_damage = sim_load.model.fields.g3d.Φ,
                                        initial_u_veloc = sim_load.model.fields.gu.u,
                                        initial_v_veloc = sim_load.model.fields.gv.v)

# Set solver params
maxiter_picard = 3
tol_picard = 1.0e-4
n_iter_viscosity = 2

solver_params = SolverParams(maxiter_picard = maxiter_picard,
                            tol_picard = tol_picard)


#parallel_spec = BasicParallelSpec()
#parallel_spec = SharedMemorySpec(ngridsx = {{run.ngridsx}}, ngridsy={{run.ngridsy}}, overlap={{run.overlap}}, niterations={{run.niterations}}) 


#Physical parameters
#we do not evolve the shelves:
evolveShelves = true

#from inversion/relaxation:
accumulation_rate = sim_load.model.fields.gh.accumulation
weertman_c = sim_load.model.params.weertman_c


params = Params(accumulation_rate = accumulation_rate,
                weertman_c = weertman_c,
                evolveShelves = evolveShelves)

 #build the model to use in the initialisation:
println("Starting to make the model")

model = Model(grid = grid,
              bed_elevation = sim_load.model.fields.gh.b, 
              params = params, 
              solver_params = solver_params,
          #    parallel_spec = parallel_spec,
              initial_conditions= initial_conditions,
#	          melt_rate = MeltRateExponentVariation(Î³T = {{ run.gammaT }}, melt_exp= {{ run.melt_exp }}, melt_partial_cell = {{ run.melt_partial_cell }}))
              melt_rate = MeltRateExponentVariation(γT = 0.0, melt_exp=2.0, melt_partial_cell = false))
     
#
println("The model is made")

 #timestepping parameters
niter0 = 0
dt = 0.005
end_time = 10.0
chkpt_freq = 25.0
pchkpt_freq = 5.0
timestepping_params = TimesteppingParams(niter0 = niter0, 
                                           dt = dt, 
                                           end_time = end_time, 
                                           chkpt_freq = chkpt_freq, 
                                           pchkpt_freq = pchkpt_freq)

##output parameters
folder = "outputs_8km_relax_test"
isdir(folder) && rm(folder, force = true, recursive = true)
mkdir(folder) #make a clean folder for outputs
outputs = (h = model.fields.gh.h,
	   b = model.fields.gh.b,
            u = model.fields.gh.u,
            v = model.fields.gh.v,
            s = model.fields.gh.s,
            h_mask = model.fields.gh.mask,
            visc_av = model.fields.gh.ηav,
            weertman_c = model.fields.gh.weertman_c,
            grounded_frac = model.fields.gh.grounded_fraction,
            av_speed = model.fields.gh.av_speed,
            melt = model.fields.gh.basal_melt,
            dhdt = model.fields.gh.dhdt,
            beta2 = model.fields.gh.β,
            accumulation = model.fields.gh.accumulation,
            bed_speed = model.fields.gh.bed_speed,
            haf = model.fields.gh.haf, 
	        uu = model.fields.gu.u,
	        vv = model.fields.gv.v,
	        gamma_T=model.melt_rate.γT,
	        melt_exp=model.melt_rate.melt_exp,
	        density_ice=model.params.density_ice,
	        density_ocean=model.params.density_ocean)
         #output velocities and thickness
             
output_freq = 1.0
output_params = OutputParams(outputs = outputs, 
                           output_path = folder,
                           output_freq = output_freq,
                           output_format = "mat",
                           zip_format = "nc")
   
                           
println("About to make simulation")
simulation = Simulation(model = model, 
                      timestepping_params = timestepping_params,
                      output_params = output_params)
           
println("The simulation is made")
#perform the simulation
run_simulation!(simulation)
println("The simulation has been run")

return model

end
