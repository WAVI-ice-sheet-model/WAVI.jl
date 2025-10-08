#This is a test of the inversion process with MISMIP_PLUS
using WAVI 
using ImageFiltering
using SparseArrays
using IterativeSolvers

function MISMIP_PLUS_inversion(simulation)

@unpack model = simulation

   #pre-select points not near data gaps:  
   outside_model_u = .!simulation.model.fields.gu.mask
   gap_mask_u = outside_model_u
   u_neargap = imfilter((gap_mask_u) .|> Int, Kernel.ones(3,3),  Pad(:replicate)) .> 0
   speed_u_mask = ((simulation.model.fields.gu.mask .== 1) .& (u_neargap .== 0))
   #mask goes to edge of domain, but mask out first and last rows.
   speed_u_mask[1, :] .= false
   speed_u_mask[end, :] .= false

   outside_model_v = .!simulation.model.fields.gv.mask
   gap_mask_v = outside_model_v
   v_neargap = imfilter((gap_mask_v) .|> Int, Kernel.ones(3,3),  Pad(:replicate)) .> 0
   speed_v_mask = ((simulation.model.fields.gv.mask .== 1) .& (v_neargap .== 0))

   data = DataFields(grid = deepcopy(simulation.model.grid),
                              speed_u = deepcopy(simulation.model.fields.gu.us),
                              speed_u_mask =  speed_u_mask,
                              speed_v = deepcopy(simulation.model.fields.gv.vs),
                              speed_v_mask =  speed_v_mask,
                              dhdt = deepcopy(simulation.model.fields.gh.dhdt),
                              accumulation_rate = deepcopy(simulation.model.params.accumulation_rate),
                              dhdtacc_mask=deepcopy(simulation.model.fields.gh.mask))

reltol=0.5
abstol=0.1
maxiter=500
gmres_restart =50
βgrounded_start=1.0e3
βfloating_start=1.0e-10
ηstart_guess = 1.0e6
βpower = 0.1
Bpower_shelf = 0.1
Bpower_grounded = 0.0
cg=false
gmres=true

inversion_params = InversionParams(reltol = reltol,
                                    maxiter = maxiter,
                                    abstol = abstol,
                                    gmres_restart = gmres_restart,
                                    βgrounded_start = βgrounded_start,
                                    βfloating_start = βfloating_start,
                                    ηstart_guess = ηstart_guess,
                                    βpower = βpower,
                                    Bpower_grounded = Bpower_grounded,
                                    Bpower_shelf = Bpower_shelf,
                                    cg = cg,
                                    gmres = gmres)
                                
#JKVstepping parameters
niter0 = 0
n_iter_out=1
max_JKV_iterations = 300
n_iter_chkpt = 1000
n_iter_pchkpt= 100

JKVstepping_params = JKVsteppingParams(niter0 = niter0, 
                                        n_iter_chkpt = n_iter_chkpt,
                                        n_iter_pchkpt = n_iter_pchkpt,
                                        n_iter_total = max_JKV_iterations,
                                        n_iter_out = n_iter_out)

inversion_output = InversionOutput()
   
dirichlet_model = deepcopy(simulation.model)

inversion = Inversion(model = dirichlet_model,
                     data = data,
                     inversion_params=inversion_params,
                     inversion_output = inversion_output)
                           
inversion_simulation = InversionSimulation(model = model, 
                                        inversion = inversion,
                                        JKVstepping_params = JKVstepping_params)
           
run_inversion_simulation!(inversion_simulation)

return inversion_simulation

end
