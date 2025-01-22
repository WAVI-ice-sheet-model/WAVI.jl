#This is a test relaxation file to see if WAVI can run  

using WAVI 
using Printf
using ImageFiltering



function inverse_driver_8km()

#
#Grid and boundary conditions
#
nx = 105
ny = 123
nσ = 12
x0 = -1825000.
y0 = -865000.
dx = 8000.0
dy = 8000.0

#h_mask=trues(nx,ny)
#u_iszero = falses(nx+1,ny); u_iszero[1,:].=true
#v_iszero=falses(nx,ny+1); v_iszero[:,1].=true; v_iszero[:,end].=true

#h_mask=Array{Float64}(undef,nx,ny);
#read!("Initial_Data/Inverse_8km_h_mask_clip.bin",h_mask)
#h_mask.=ntoh.(h_mask)

h_mask=Array{Float64}(undef,nx,ny);
read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_github/WAVI-WAIS-setups/inversion_data/bedmachinev3/full_stripe_fix_8km/Inverse_8km_h_mask_clip_BedmachineV3_FULL_stripe_fix.bin",h_mask)
h_mask.=ntoh.(h_mask)

u_iszero=Array{Float64}(undef,nx+1,ny);
read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_github/WAVI-WAIS-setups/inversion_data/bedmachinev3/full_stripe_fix_8km/Inverse_8km_uiszero_clip_BedmachineV3_FULL_stripe_fix.bin",u_iszero)
u_iszero.=ntoh.(u_iszero)

v_iszero=Array{Float64}(undef,nx,ny+1);
read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_github/WAVI-WAIS-setups/inversion_data/bedmachinev3/full_stripe_fix_8km/Inverse_8km_viszero_clip_BedmachineV3_FULL_stripe_fix.bin",v_iszero)
v_iszero.=ntoh.(v_iszero)

sigma_grid=Array{Float64}(undef,nσ);
read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_github/WAVI-WAIS-setups/inversion_data/bedmachinev3/full_stripe_fix_8km/Inverse_8km_sigma_grid_BedmachineV3_FULL_stripe_fix.bin",sigma_grid)
sigma_grid.=ntoh.(sigma_grid)

basin_ID=Array{Float64}(undef,nx,ny);
read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_github/WAVI-WAIS-setups/inversion_data/bedmachinev3/full_stripe_fix_8km/Inverse_8km_basinID_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",basin_ID)
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
read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_github/WAVI-WAIS-setups/inversion_data/bedmachinev3/full_stripe_fix_8km/Inverse_8km_bed_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",bed)
bed.=ntoh.(bed)

h=Array{Float64}(undef,nx,ny);
read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_github/WAVI-WAIS-setups/inversion_data/bedmachinev3/full_stripe_fix_8km/Inverse_8km_thickness_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",h)
h.=ntoh.(h) 

 viscosity=Array{Float64}(undef,nx,ny,nσ);
 read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_github/WAVI-WAIS-setups/inversion_data/bedmachinev3/full_stripe_fix_8km/Inverse_8km_viscosity3D_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",viscosity)
 viscosity.=ntoh.(viscosity)

 temp=Array{Float64}(undef,nx,ny,nσ);
 read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_github/WAVI-WAIS-setups/inversion_data/bedmachinev3/full_stripe_fix_8km/Inverse_8km_3Dtemp_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",temp)
 temp.=ntoh.(temp)

 damage=Array{Float64}(undef,nx,ny,nσ);
  read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_github/WAVI-WAIS-setups/inversion_data/bedmachinev3/full_stripe_fix_8km/Inverse_8km_damage3D_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",damage)
 damage.=ntoh.(damage)

 
 weertman_c=Array{Float64}(undef,nx,ny);
 read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_github/WAVI-WAIS-setups/inversion_data/bedmachinev3/full_stripe_fix_8km/Inverse_8km_WeertmanC_clip_adjusted_noNan_BedmachineV3_FULL_stripe_fix.bin",weertman_c)
 weertman_c.=ntoh.(weertman_c)

 accumulation_rate=Array{Float64}(undef,nx,ny);
 read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_github/WAVI-WAIS-setups/inversion_data/bedmachinev3/full_stripe_fix_8km/Inverse_8km_accumulation_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",accumulation_rate)
 accumulation_rate.=ntoh.(accumulation_rate)


 dhdt=Array{Float64}(undef,nx,ny);
 read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_github/WAVI-WAIS-setups/inversion_data/bedmachinev3/full_stripe_fix_8km/Inverse_8km_dhdt_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",dhdt)
 dhdt.=ntoh.(dhdt)

 gu_u=Array{Float64}(undef,nx+1,ny);
 read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_github/WAVI-WAIS-setups/inversion_data/bedmachinev3/full_stripe_fix_8km/Inverse_8km_u_velocs_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",gu_u)
 gu_u.=ntoh.(gu_u)

 gv_v=Array{Float64}(undef,nx,ny+1);
 read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_github/WAVI-WAIS-setups/inversion_data/bedmachinev3/full_stripe_fix_8km/Inverse_8km_v_velocs_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",gv_v)
 gv_v.=ntoh.(gv_v)

initial_conditions = InitialConditions(initial_thickness = h,
                                        initial_viscosity = viscosity,
                                        initial_temperature = temp,
                                        initial_damage = damage,
                                        initial_u_veloc = gu_u,
                                        initial_v_veloc = gv_v)


#Can I add in initialparams here? So that I can vary these with the ensembler in the end??

#solver parameters
#
maxiter_picard = 10
tol_picard = 1.0e-4
solver_params = SolverParams(maxiter_picard = maxiter_picard,
                            tol_picard = tol_picard)

#
#Physical parameters
#default_thickness = 100.0 #set the initial condition this way
#accumulation_rate = 0.3
#we do not evolve the shelves:
evolveShelves = false

accumulation_rate = accumulation_rate



params = Params(accumulation_rate = accumulation_rate,
                weertman_c = weertman_c,
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

#timestepping parameters
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


@printf "The inversion test is about to be done"
#perform the inversion

#First read in the data to be used for the inversion:

accumulation_rate=Array{Float64}(undef,nx,ny);
read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_github/WAVI-WAIS-setups/inversion_data/bedmachinev3/full_stripe_fix_8km/Inverse_8km_accumulation_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",accumulation_rate)
accumulation_rate.=ntoh.(accumulation_rate)

dhdt=Array{Float64}(undef,nx,ny);
read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_github/WAVI-WAIS-setups/inversion_data/bedmachinev3/full_stripe_fix_8km/Inverse_8km_dhdt_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",dhdt)
dhdt.=ntoh.(dhdt)

dhdtaccmask=Array{Float64}(undef,nx,ny);
read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_github/WAVI-WAIS-setups/inversion_data/bedmachinev3/full_stripe_fix_8km/Inverse_8km_dhdtaccmask_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",dhdtaccmask)
dhdtaccmask.=ntoh.(dhdtaccmask)

udata=Array{Float64}(undef,nx+1,ny);
read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_github/WAVI-WAIS-setups/inversion_data/bedmachinev3/full_stripe_fix_8km/Inverse_8km_udata_clip_BedmachineV3_FULL_stripe_fix.bin",udata)
udata.=ntoh.(udata)

vdata=Array{Float64}(undef,nx,ny+1);
read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_github/WAVI-WAIS-setups/inversion_data/bedmachinev3/full_stripe_fix_8km/Inverse_8km_vdata_clip_BedmachineV3_FULL_stripe_fix.bin",vdata)
vdata.=ntoh.(vdata)

udatamask=Array{Float64}(undef,nx+1,ny);
read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_github/WAVI-WAIS-setups/inversion_data/bedmachinev3/full_stripe_fix_8km/Inverse_8km_udatamask_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",udatamask)
udatamask.=ntoh.(udatamask)

vdatamask=Array{Float64}(undef,nx,ny+1);
read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_github/WAVI-WAIS-setups/inversion_data/bedmachinev3/full_stripe_fix_8km/Inverse_8km_vdatamask_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",vdatamask)
vdatamask.=ntoh.(vdatamask)

udatamask_combo = ((udatamask .== 1) .& (model.fields.gu.mask .== 1))
vdatamask_combo = ((vdatamask .== 1) .& (model.fields.gv.mask .== 1))
dhdtaccmask_combo = ((dhdtaccmask .== 1) .& (model.fields.gh.mask .== 1))

udatamask_combo = convert(Array{Bool,2}, udatamask_combo)
vdatamask_combo = convert(Array{Bool,2}, vdatamask_combo)
dhdtaccmask_combo = convert(Array{Bool,2},dhdtaccmask_combo)

inversion = Inversion(grid = grid,
                    speed_u = udata,
                    speed_u_mask = udatamask_combo,
                    speed_v = vdata,
                    speed_v_mask = vdatamask_combo,
                    dhdt = dhdt,
                    accumulation_rate = accumulation_rate,
                    dhdtacc_mask=dhdtaccmask_combo)

#Then run the inversion:
run_inversion!(simulation,inversion)

@printf "The inversion test has been run"


return model


end


#= 
"""
    get_u_mask(h_mask)

Find mask of valid grid points on u-grid corresponding to a mask defined on h-grid.

"""
function get_u_mask(h_mask)
    #include all u faces next to a selected center
    (nx,ny)=size(h_mask)
    u_mask=falses(nx+1,ny)
    u_mask[1:end-1,1:end]=u_mask[1:end-1,1:end].|h_mask
    u_mask[2:end,1:end]=u_mask[2:end,1:end].|h_mask
    return u_mask
end

 u_mask = get_u_mask(h_mask)

    #u-grid
    gu=UGrid(
    nxu=grid.nx+1,
    nyu=grid.ny,
    dx=grid.dx,
    dy=grid.dy,
    mask=u_mask,
    u_isfixed=grid.u_isfixed,
    u=deepcopy(initial_conditions.initial_u_veloc),
    levels=solver_params.levels
    ) =#