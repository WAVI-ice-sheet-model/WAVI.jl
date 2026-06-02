#This is a test 8km relaxation script: 

using WAVI 

function relax_driver_8km()

#Load up the model from the inversion:
n_iter_inversion=5;
  
file_path="/data/glacio/icesheets/shelf_thin/wip/chll1/WAVI_inversion/WAVI.jl/outputs_8km_inversion_test_speed_store"
#PChkpt_0000000015.jld2"

    n_iter_string =  lpad(n_iter_inversion, 10, "0"); #filename as a string with 10 digits
        println("Requested inversion is iteration $n_iter_inversion. Looking for inversion file...")
    #    try 
            #for running on vscode vis workstation:
            filename = string(file_path, "/PChkpt_", n_iter_string, ".jld2")
            #for running with ensembler:
            #filename = string("PChkpt_", n_iter_string, ".jld2")
            println("Looking for file: ", filename)

            if isfile(filename)
                println("inversion file exists")
                sim_load = load(filename, "inversion_simulation")
                println("Imversion pickup successful")
            else
                error("File $filename not found!")
            end
            

#load up grid used for inversion:  
grid = sim_load.model.grid
          
# set up initial conditions from inversion:        
initial_conditions = InitialConditions(initial_thickness =  sim_load.model.fields.gh.h,
                                        initial_viscosity = sim_load.model.fields.g3d.η,
                                        initial_temperature = sim_load.model.fields.g3d.θ,
                                        initial_damage = sim_load.model.fields.g3d.Φ)
#                                        initial_u_veloc = gu_u,
 #                                       initial_v_veloc = gv_v)

# Set solver params
maxiter_picard = 1
tol_picard = 1.0e-4
n_iter_viscosity = 2

solver_params = SolverParams(maxiter_picard = maxiter_picard,
                            tol_picard = tol_picard,
                            n_iter_viscosity = n_iter_viscosity)


#Physical parameters
#we do not evolve the shelves:
evolveShelves = false

#from inversion:
accumulation_rate = sim_load.model.fields.gh.accumulation - sim_load.inversion.data_fields.ghdata.dhdt

#calculate weertman_c directly from inversion:
weertman_c_newly_grounded=1.0e4
weertman_c=zeros(size(sim_load.model.fields.gh.β))
weertman_c[sim_load.model.fields.gh.mask] = sim_load.model.fields.gh.β[sim_load.model.fields.gh.mask] ./(( sqrt.(sim_load.model.fields.gh.bed_speed[sim_load.model.fields.gh.mask].^2 .+  sim_load.model.params.weertman_reg_speed^2)).^(1.0./sim_load.model.params.weertman_m .- 1.0))
#weertman_c = sim_load.model.fields.gh.β ./( sqrt.(sim_load.model.fields.gh.bed_speed.^2 .+  sim_load.model.params.weertman_reg_speed^2)).^(1.0./sim_load.model.params.weertman_m .- 1.0)

#adjust for small value at partially grounded cells:
weertman_c_adjusted=zeros(size(sim_load.model.fields.gh.β))
aground=(sim_load.model.fields.gh.haf .>=0)

f1 = findall(aground .== 1)
f2 = findall(aground .== 0)

weertman_c_adjusted[f1] .= weertman_c[f1]./sim_load.model.fields.gh.grounded_fraction[f1]
weertman_c_adjusted[f2] .= weertman_c_newly_grounded
# remove any nans
weertman_c_adjusted[isnan.(weertman_c_adjusted) .& .!(weertman_c_adjusted .∈ Ref(sim_load.model.fields.gh.mask))] .= -9999.0

params = Params(accumulation_rate = accumulation_rate,
                weertman_c = weertman_c_adjusted,
                evolveShelves = evolveShelves)


 #build the model to use in the initialisation:
println("Starting to make the model")
model = Model(grid = grid,
              bed_elevation = sim_load.model.fields.gh.b, 
              params = params, 
              solver_params = solver_params,
              initial_conditions= initial_conditions) 
#
println("The model is made")

 #timestepping parameters
niter0 = 0
dt = 0.1
end_time = 2.0
chkpt_freq = 25.0
pchkpt_freq = 1.0
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
            u = model.fields.gh.u,
            v = model.fields.gh.v,
            b = model.fields.gh.b,
            s = model.fields.gh.s,
            h_mask = model.fields.gh.mask)
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
