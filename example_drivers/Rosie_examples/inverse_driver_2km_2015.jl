#This is a test inversion file to see if WAVI can run  

using WAVI 
using Printf
using ImageFiltering

function inverse_driver_2km_2015()

#
#Grid and boundary conditions
#
nx = 398
ny = 468
nσ = 26
x0 = -1789000.
y0 = -833000.
dx = 2000.0
dy = 2000.0

h_mask=Array{Float64}(undef,nx,ny);
read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_2015/Mask_0_2km/h_mask.bin",h_mask)
h_mask.=ntoh.(h_mask)

u_iszero=Array{Float64}(undef,nx+1,ny);
read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_2015/Mask_0_2km/u_iszero.bin",u_iszero)
u_iszero.=ntoh.(u_iszero)

v_iszero=Array{Float64}(undef,nx,ny+1);
read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_2015/Mask_0_2km/v_iszero.bin",v_iszero)
v_iszero.=ntoh.(v_iszero)

sigma_grid=Array{Float64}(undef,nσ);
read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_2015/Mask_0_2km/sigma_grid.bin",sigma_grid)
sigma_grid.=ntoh.(sigma_grid)

basin_ID=Array{Float64}(undef,nx,ny);
read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_2015/Mask_0_2km/basinID.bin",basin_ID)
basin_ID.=ntoh.(basin_ID)


grid = Grid(nx = nx, 
            ny = ny,   
            nσ = nσ, 
            x0 = x0, 
            y0 = y0, 
            dx = dx, 
            dy = dy,
            h_mask = h_mask, 
            u_iszero = u_iszero, 
            v_iszero = v_iszero,
            σ = sigma_grid,
            basin_ID = basin_ID)
#
#Bed 
#
bed=Array{Float64}(undef,nx,ny);
read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_2015/Mask_0_2km/bed.bin",bed)
bed.=ntoh.(bed)

h=Array{Float64}(undef,nx,ny);
read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_2015/Mask_0_2km/thickness.bin",h)
h.=ntoh.(h) 

#this should read in from inversion, because is necessary for julia setup,but is it used???
# viscosity=Array{Float64}(undef,nx,ny,nσ);
 #read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_2015/Mask_0_2km/Inverse_2km_viscosity3D_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",viscosity)
 #viscosity.=ntoh.(viscosity)

 temp=Array{Float64}(undef,nx,ny,nσ);
 read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_2015/Mask_0_2km/temps.bin",temp)
 temp.=ntoh.(temp)

 #check this comes in as all zeros...
# damage=Array{Float64}(undef,nx,ny,nσ);
 # read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_2015/Mask_0_2km/Inverse_2km_damage3D_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",damage)
# damage.=ntoh.(damage)

#all zeros and isn't used?? Does it actually need to be here?
 #weertman_c=Array{Float64}(undef,nx,ny);
 #read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_2015/Mask_0_2km/Inverse_2km_WeertmanC_clip_adjusted_noNan_BedmachineV3_FULL_stripe_fix.bin",weertman_c)
 #weertman_c.=ntoh.(weertman_c)

 #same as accumulation data? Redundant??
 accumulation_rate=Array{Float64}(undef,nx,ny);
 read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_2015/Mask_0_2km/accumulation_data.bin",accumulation_rate)
 accumulation_rate.=ntoh.(accumulation_rate)

#This is from data
 dhdt=Array{Float64}(undef,nx,ny);
 read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_2015/Mask_0_2km/dhdt_data.bin",dhdt)
 dhdt.=ntoh.(dhdt)

 #Hope these aren't actually needed...
 #gu_u=Array{Float64}(undef,nx+1,ny);
 #read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_2015/Mask_0_2km/Inverse_2km_u_velocs_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",gu_u)
 #gu_u.=ntoh.(gu_u)
#
 #gv_v=Array{Float64}(undef,nx,ny+1);
 #read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_2015/Mask_0_2km/Inverse_2km_v_velocs_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",gv_v)
 #gv_v.=ntoh.(gv_v)

initial_conditions = InitialConditions(initial_thickness = h,
#                                        initial_viscosity = viscosity,
                                        initial_temperature = temp)
 #                                       initial_damage = damage)
#                                        initial_u_veloc = gu_u,
 #                                       initial_v_veloc = gv_v)

#solver parameters
#
maxiter_picard = 1
tol_picard = 1.0e-4
n_iter_viscosity = 2

solver_params = SolverParams(maxiter_picard = maxiter_picard,
                            tol_picard = tol_picard,
                            n_iter_viscosity = n_iter_viscosity)


#
#Physical parameters
#default_thickness = 100.0 #set the initial condition this way
#accumulation_rate = 0.3
#we do not evolve the shelves:
evolveShelves = false

accumulation_rate = accumulation_rate

params = Params(accumulation_rate = accumulation_rate,
                weertman_c = 0.0,
                evolveShelves = evolveShelves)


 #build the model to use in the initialisation:
@printf "Starting to make the model"
model = Model(grid = grid,
              bed_elevation = bed, 
              params = params, 
              solver_params = solver_params,
              initial_conditions= initial_conditions) 
#
@sprintf "The model is made"

#But now I don't want to do any timestepping, I want to initialise instead and produce new model, no simulation...

#= #timestepping parameters
niter0 = 0
dt = 0.1
end_time = 1.0
chkpt_freq = 25.0
pchkpt_freq = 25.0
timestepping_params = TimesteppingParams(niter0 = niter0, 
                                           dt = dt, 
                                           end_time = end_time, 
                                           chkpt_freq = chkpt_freq, 
                                           pchkpt_freq = pchkpt_freq)

##output parameters
folder = "outputs_2km_relax_test"
isdir(folder) && rm(folder, force = true, recursive = true)
mkdir(folder) #make a clean folder for outputs
outputs = (h = model.fields.gh.h,
            u = model.fields.gh.u,
            v = model.fields.gh.v,
            b = model.fields.gh.b,
            s = model.fields.gh.s,
            h_mask = model.fields.gh.mask)
             #output velocities and thickness
output_freq = 0.1
output_params = OutputParams(outputs = outputs, 
                           output_path = folder,
                           output_freq = output_freq,
                           output_format = "mat",
                           zip_format = "nc")
   
                           
@printf "About to make simulation"
simulation = Simulation(model = model, 
                      timestepping_params = timestepping_params,
                      output_params = output_params)
           
 @printf "The simulation is made"
#perform the simulation
#run_simulation!(simulation)
@printf "The simulation has NOT been run" 
 =#

@printf "The inversion test is about to be done"
#perform the inversion

#First read in the data to be used for the inversion:

#= accumulation_rate=Array{Float64}(undef,nx,ny);
read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_2015/Mask_0_2km/Inverse_2km_accumulation_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",accumulation_rate)
accumulation_rate.=ntoh.(accumulation_rate)

dhdt=Array{Float64}(undef,nx,ny);
read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_2015/Mask_0_2km/Inverse_2km_dhdt_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",dhdt)
dhdt.=ntoh.(dhdt) =#

dhdtaccmask=Array{Float64}(undef,nx,ny);
read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_2015/Mask_0_2km/dhdt_acc_mask.bin",dhdtaccmask)
dhdtaccmask.=ntoh.(dhdtaccmask)

udata=Array{Float64}(undef,nx+1,ny);
read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_2015/Mask_0_2km/udata.bin",udata)
udata.=ntoh.(udata)

vdata=Array{Float64}(undef,nx,ny+1);
read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_2015/Mask_0_2km/vdata.bin",vdata)
vdata.=ntoh.(vdata)

udatamask=Array{Float64}(undef,nx+1,ny);
read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_2015/Mask_0_2km/udata_mask.bin",udatamask)
udatamask.=ntoh.(udatamask)

vdatamask=Array{Float64}(undef,nx,ny+1);
read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_2015/Mask_0_2km/vdata_mask.bin",vdatamask)
vdatamask.=ntoh.(vdatamask)

#pre-select points not near data gaps:  
u_neargap = imfilter((model.fields.gu.mask .& .!(udatamask .> 0)) .|> Int, Kernel.ones(3,3), Pad(1,1)) .> 0
v_neargap = imfilter((model.fields.gv.mask .& .!(vdatamask .> 0)) .|> Int, Kernel.ones(3,3), Pad(1,1)) .> 0
#data_mask_u = .!gudata.neargap .& .!gu.u_isfixed .& gu.mask

    
udatamask_combo = ((udatamask .== 1) .& (model.fields.gu.mask .== 1) .& (u_neargap .== 0))
vdatamask_combo = ((vdatamask .== 1) .& (model.fields.gv.mask .== 1) .& (v_neargap .== 0))
dhdtaccmask_combo = ((dhdtaccmask .== 1) .& (model.fields.gh.mask .== 1))
#
udatamask_combo = convert(Array{Bool,2}, udatamask_combo)
vdatamask_combo = convert(Array{Bool,2}, vdatamask_combo)
dhdtaccmask_combo = convert(Array{Bool,2},dhdtaccmask_combo)

println("nnz inmodel.fields.gu.mask is " ,count(!iszero, model.fields.gu.mask))
println("nnz in udatamask is " ,count(!iszero, udatamask))
println("nnz in udatamask_combo is " ,count(!iszero, udatamask_combo))
println("nnz in vdatamask_combo is " ,count(!iszero, vdatamask_combo))

gmres_reltol=0.5
gmres_abstol=0.01
gmres_maxiter=500
gmres_restart =100
βgrounded_start=1.e4
βfloating_start=1.0e-10
ηstart_guess = 1.0e7
βpower = 0.1
Bpower_shelf = 0.1
Bpower_grounded = 0.01
#max_JKV_iterations = 5

inversion_params = InversionParams(gmres_reltol = gmres_reltol,
                                    gmres_maxiter = gmres_maxiter,
                                    gmres_restart = gmres_restart,
                                    βgrounded_start = βgrounded_start,
                                    βfloating_start = βfloating_start,
                                    ηstart_guess = ηstart_guess,
                                    βpower = βpower)
               #                     max_JKV_iterations = max_JKV_iterations)

                                
#JKVstepping parameters
niter0 = 0
n_iter_out=1
max_JKV_iterations = 30
n_iter_chkpt = 100
n_iter_pchkpt= 5

JKVstepping_params = JKVsteppingParams(niter0 = niter0, 
                                        n_iter_chkpt = n_iter_chkpt,
                                        n_iter_pchkpt = n_iter_pchkpt,
                                        n_iter_total = max_JKV_iterations,
                                        n_iter_out = n_iter_out)

JKV=zeros(max_JKV_iterations)
JRMS=zeros(max_JKV_iterations)

inversion_output = InversionOutput(JKV=JKV,
                                    JRMS=JRMS)

                
inversion = Inversion(grid = grid,
                    bed_elevation=bed,
                    inversion_params=inversion_params,
                    solver_params=solver_params,
                    speed_u = udata,
                    speed_u_mask = udatamask_combo,
                    speed_v = vdata,
                    speed_v_mask = vdatamask_combo,
                    dhdt = dhdt,
                    accumulation_rate = accumulation_rate,
                    dhdtacc_mask=dhdtaccmask_combo,
                    initial_conditions=initial_conditions,
                    params = params,
                    inversion_output=inversion_output)
 #                   melt_rate=nothing,
 #                   parallel_spec=nothing)

 @printf "About to make inversion_simulation"

 ##output parameters
folder = "Y2015_Mask_0_2km_inversion_slim_more_restarts"
isdir(folder) && rm(folder, force = true, recursive = true)
mkdir(folder) #make a clean folder for outputs
outputs = (h = model.fields.gh.h,
            u = model.fields.gu.u,
            uh = model.fields.gh.u,
            v = model.fields.gv.v,
            vh = model.fields.gh.v,
            b = model.fields.gh.b,
            s = model.fields.gh.s,
            av_speed = model.fields.gh.av_speed,
            av_speed_d = inversion.fields.gh.av_speed,
            dhdt_n = model.fields.gh.dhdt,
            h_mask = model.fields.gh.mask,
            u_mask =  model.fields.gu.mask,
            v_mask =  model.fields.gv.mask,
            grounded_fraction = model.fields.gh.grounded_fraction,
            beta = model.fields.gh.β,
            beta_eff = model.fields.gh.βeff,
            preBfactor = model.fields.gh.preBfactor,
            quad_f0 = model.fields.gh.quad_f0,
            quad_f1 = model.fields.gh.quad_f1,
            quad_f2 = model.fields.gh.quad_f2,
            eta_av = model.fields.gh.ηav,
            udata = inversion.data_fields.gudata.speed_u,
            vdata = inversion.data_fields.gvdata.speed_v,
            udata_mask = inversion.data_fields.gudata.mask,
            vdata_mask = inversion.data_fields.gvdata.mask,
            dhdt_data = inversion.data_fields.ghdata.dhdt,
            dhdtaccdata_mask = inversion.data_fields.ghdata.mask,
            u_d = inversion.fields.gu.u,
            u_dh = inversion.fields.gh.u,
            v_d = inversion.fields.gv.v,
            v_dh = inversion.fields.gh.v,
            u_sd = inversion.fields.gu.us,
            v_sd = inversion.fields.gv.vs,
            u_s = model.fields.gu.us,
            v_s = model.fields.gv.vs,
            u_sd_h = inversion.fields.gh.us,
            v_sd_h = inversion.fields.gh.vs,
            u_s_h = model.fields.gh.us,
            v_s_h = model.fields.gh.vs,
            u_bd_h = inversion.fields.gh.ub,
            v_bd_h = inversion.fields.gh.vb,
            u_b_h = model.fields.gh.ub,
            v_b_h = model.fields.gh.vb,
            tau_surf_u= inversion.fields.gu.τsurf,
            tau_surf_v = inversion.fields.gv.τsurf,
            sigmazz_h= inversion.fields.gh.σzzsurf,
            JKV=inversion.inversion_output.JKV,
            JRMS=inversion.inversion_output.JRMS,
            ud_resid=inversion.fields.gu.residual,
            vd_resid=inversion.fields.gv.residual,
            udata_resid=inversion.data_fields.gudata.residual,
            vdata_resid=inversion.data_fields.gvdata.residual,   
            hdata_resid=inversion.data_fields.ghdata.residual,
            shelf_strain_rate=model.fields.gh.shelf_strain_rate,
            shelf_strain_rate_d=inversion.fields.gh.shelf_strain_rate
           # gh=model.fields.gh,
            )

output_freq = 1
output_params = OutputParams(outputs = outputs, 
                           output_path = folder,
                           output_freq = output_freq,
                           output_format = "mat",
                           zip_format = "nc")
   

inversion_simulation = InversionSimulation(model = model, 
                                        inversion = inversion,
                                        JKVstepping_params = JKVstepping_params,
                                        output_params = output_params)
           

#Then run the inversion:
run_inversion_simulation!(inversion_simulation)

#run_inversion!(simulation,inversion)

@printf "The inversion test has been run"

return model

end
