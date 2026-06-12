
using WAVI, NCDatasets, Test, Parameters, Interpolations
using WAVI.Parameters


"""
           iceberg_grid(; nx, ny)

       Build the small Iceberg test grid with fixed velocities on the centre lines.
"""
function iceberg_grid(; nx = 12, ny = 8)
    u_isfixed = falses(nx + 1, ny)
    u_isfixed[div(nx, 2) + 1, :] .= true
    v_isfixed = falses(nx, ny + 1)
    v_isfixed[:, div(ny, 2) + 1] .= true
    return Grid(
        nx = nx,
        ny = ny,
        nσ = 4,
        x0 = -50000.0,
        y0 = -50000.0,
        dx = 5000.0,
        dy = 5000.0,
        u_isfixed = u_isfixed,
        v_isfixed = v_isfixed,
           )
end

function climate_forcing_test(grid)

params = Params()

#random thickness iceberg, but not too thin
thickness_scale = 500.0
reference_thickness = thickness_scale .* (0.1 .+ 0.9*rand(grid.nx,grid.ny))

nzforcing = 3
z_forcing = .- thickness_scale * 0.9 .* [0.0,0.5,1.0]     

#deep enough to ensure floatation
bed_elevation = -1000*ones(grid.nx,grid.ny)

ISMIP7_config = ISMIP7(gcm = "TEST",scenario = "MOCK")

shelf_melt_rate = ISMIP7MeltRate(;ISMIP7_config,z_forcing)

#Just to get default parameters for later
@unpack K, shelf_slope, ρ_ocean, ρ_ice, c_ocean, L_ice, β_s, g, f, melt_partial_cell = shelf_melt_rate 
C_melt = params.sec_per_year*K*shelf_slope*(ρ_ocean/ρ_ice)*g*β_s*(c_ocean/L_ice)^2/(2*abs(f)) 

manufactured_melt_solution = 0.01*ones(grid.nx,grid.ny,nzforcing,2)

delta = (one(typeof(ρ_ice/ρ_ocean)) .- ρ_ice./ρ_ocean)

reference_elevation = reference_thickness * delta  
reference_smb = rand(grid.nx,grid.ny)

mask = Int8.(rand(Bool,grid.nx,grid.ny,2))
NCDataset("ice_shelf_collapse_mask_TEST_sspMOCK_ismip7_5000m_0.nc","c") do ds
    defVar(ds,"mask",mask[:,:,1],("x","y","t"))
end
NCDataset("ice_shelf_collapse_mask_TEST_sspMOCK_ismip7_5000m_1.nc","c") do ds
    defVar(ds,"mask",mask[:,:,2],("x","y","t"))
end
acabf_anomaly = rand(Float32,grid.nx,grid.ny,2)
NCDataset("acabf-anomaly_AIS_TEST_sspMOCK_SDBN1-5000m_v2_0_yearlyaveraged.nc","c") do ds
    defVar(ds,"acabf-anomaly",acabf_anomaly[:,:,1],("x","y","t"))
end
NCDataset("acabf-anomaly_AIS_TEST_sspMOCK_SDBN1-5000m_v2_1_yearlyaveraged.nc","c") do ds
    defVar(ds,"acabf-anomaly",acabf_anomaly[:,:,2],("x","y","t"))
end
dacabfdz = rand(Float32,grid.nx,grid.ny,2)
NCDataset("dacabfdz_AIS_TEST_sspMOCK_SDBN1-5000m_v2_0_yearlyaveraged.nc","c") do ds
    defVar(ds,"dacabfdz",dacabfdz[:,:,1],("x","y","t"))
end
NCDataset("dacabfdz_AIS_TEST_sspMOCK_SDBN1-5000m_v2_1_yearlyaveraged.nc","c") do ds
    defVar(ds,"dacabfdz",dacabfdz[:,:,2],("x","y","t"))
end

salinity_floor = 0.1
so = rand(Float32,grid.nx,grid.ny,nzforcing,2) .+ salinity_floor
NCDataset("so_AIS_TEST_sspMOCK_ocean_v3_5000m_0.nc","c") do ds
    defVar(ds,"so",so[:,:,:,1],("x","y","z","t"))
end

NCDataset("so_AIS_TEST_sspMOCK_ocean_v3_5000m_1.nc","c") do ds
    defVar(ds,"so",so[:,:,:,2],("x","y","z","t"))
end

thetao = rand(Float32,grid.nx,grid.ny,nzforcing,2)
NCDataset("thetao_AIS_TEST_sspMOCK_ocean_v3_5000m_0.nc","c") do ds
    defVar(ds,"thetao",thetao[:,:,:,1],("x","y","z","t"))
end
NCDataset("thetao_AIS_TEST_sspMOCK_ocean_v3_5000m_1.nc","c") do ds
    defVar(ds,"thetao",thetao[:,:,:,2],("x","y","z","t"))
end

#tf = rand(Float32,grid.nx,grid.ny,nzforcing,2)
#Choose tf to match manufactured melt solution.
tf = abs.(manufactured_melt_solution./(C_melt*so)).^-0.5.*(manufactured_melt_solution./(C_melt*so))
NCDataset("tf_AIS_TEST_sspMOCK_ocean_v3_5000m_0.nc","c") do ds
    defVar(ds,"tf",tf[:,:,:,1],("x","y","z","t"))
end
NCDataset("tf_AIS_TEST_sspMOCK_ocean_v3_5000m_1.nc","c") do ds
    defVar(ds,"tf",tf[:,:,:,2],("x","y","z","t"))
end



surface_mass_balance = ISMIP7SMB(;ISMIP7_config,reference_elevation,reference_smb)
fracture = ISMIP7Hydrofracture(;ISMIP7_config)

initial_conditions = InitialConditions(initial_thickness = reference_thickness)

model = Model(grid,bed_elevation;params,initial_conditions,shelf_melt_rate,fracture,surface_mass_balance)

timestepping_params=TimesteppingParams(dt=0.1,n_iter_total=20,ntimesteps_climate_forcing_update=10)

    outputs = (
        h = "fields.gh.h",
        s = "fields.gh.s",
        u = "fields.gh.u",
        v = "fields.gh.v",
        b = "fields.gh.b",
        basal_melt = "fields.gh.basal_melt",
        accumulation = "fields.gh.accumulation",
        reference_elevation = "surface_mass_balance.reference_elevation",
        vertical_smb_gradient = "surface_mass_balance.vertical_smb_gradient",
        smb_anomaly = "surface_mass_balance.smb_anomaly",
        reference_smb = "surface_mass_balance.reference_smb",
        ice_shelf_collapse_mask = "fracture.ice_shelf_collapse_mask",
        damage = "fields.g3d.Φ",
        so = "shelf_melt_rate.S_loc",
        thetao = "shelf_melt_rate.T_loc",
        tf = "shelf_melt_rate.Tf_loc"
    )

output_freq = 1

zip_format = "nc"

output_params = OutputParams(outputs;output_freq,zip_format = zip_format)

simulation = Simulation(;model,timestepping_params,output_params)

run_simulation!(simulation)


ds = NCDataset("outfile.nc")

damage1 = load("outfile0000000010.jld2")["damage"]
damage2 = load("outfile0000000020.jld2")["damage"]


@test ds["reference_elevation"][:,:,1] == reference_elevation
@test ds["reference_elevation"][:,:,2] == reference_elevation

@test ds["reference_smb"][:,:,1] == reference_smb
@test ds["reference_smb"][:,:,2] == reference_smb

@test ds["smb_anomaly"][:,:,1] == acabf_anomaly[:,:,1]
@test ds["smb_anomaly"][:,:,2] == acabf_anomaly[:,:,2]

@test ds["accumulation"][:,:,1] ≈ ((ds["s"][:,:,1].- reference_elevation).*dacabfdz[:,:,1] .+ acabf_anomaly[:,:,1] .+ reference_smb)
@test ds["accumulation"][:,:,2] ≈ ((ds["s"][:,:,2].- reference_elevation).*dacabfdz[:,:,2] .+ acabf_anomaly[:,:,2] .+ reference_smb)

@test ds["ice_shelf_collapse_mask"][:,:,1] == mask[:,:,1]
@test ds["ice_shelf_collapse_mask"][:,:,2] == mask[:,:,2]

@test all(damage1[findall(Bool.(ds["ice_shelf_collapse_mask"][:,:,1] .!= 0.0)),:] .== fracture.damage_value)
@test all(damage2[findall(Bool.(ds["ice_shelf_collapse_mask"][:,:,2] .!= 0.0)),:] .== fracture.damage_value)

so_from_model1 = load("outfile0000000010.jld2")["so"]
so_from_model2 = load("outfile0000000020.jld2")["so"]
@test so_from_model1 == so[:,:,:,1]
@test so_from_model2 == so[:,:,:,2]

thetao_from_model1 = load("outfile0000000010.jld2")["thetao"]
thetao_from_model2 = load("outfile0000000020.jld2")["thetao"]
@test thetao_from_model1 == thetao[:,:,:,1]
@test thetao_from_model2 == thetao[:,:,:,2]

tf_from_model1 = load("outfile0000000010.jld2")["tf"]
tf_from_model2 = load("outfile0000000020.jld2")["tf"]
@test tf_from_model1 == tf[:,:,:,1]
@test tf_from_model2 == tf[:,:,:,2]

basal_melt1=ds["basal_melt"][:,:,1] 
basal_melt2=ds["basal_melt"][:,:,2] 

expected_melt_rate_3d_time1 = zeros(grid.nx,grid.ny,nzforcing)
expected_melt_rate_3d_time1 .= C_melt.*so[:,:,:,1].*tf[:,:,:,1].*abs.(tf[:,:,:,1])
expected_melt_rate_3d_time2 = zeros(grid.nx,grid.ny,nzforcing)
expected_melt_rate_3d_time2 .= C_melt.*so[:,:,:,2].*tf[:,:,:,2].*abs.(tf[:,:,:,2])

zb1 = - (ρ_ice/ρ_ocean) * ds["h"][:,:,1]
zb2 = - (ρ_ice/ρ_ocean) * ds["h"][:,:,2]

iarr = [i for i = 1:grid.nx]
jarr = [j for j = 1:grid.ny]
itp1 = interpolate((iarr,jarr,z_forcing[end:-1:1]),expected_melt_rate_3d_time1[:,:,end:-1:1],(NoInterp(),NoInterp(),Gridded(Constant())))
itp1 = extrapolate(itp1,Flat())
expected_melt_rate_2d_time1 = itp1.([i for i = 1:grid.nx,j=1:grid.ny],[j for i = 1:grid.nx,j=1:grid.ny],zb1)

itp2 = interpolate((iarr,jarr,z_forcing[end:-1:1]),expected_melt_rate_3d_time2[:,:,end:-1:1],(NoInterp(),NoInterp(),Gridded(Constant())))
itp2 = extrapolate(itp2,Flat())
expected_melt_rate_2d_time2 = itp2.([i for i = 1:grid.nx,j=1:grid.ny],[j for i = 1:grid.nx,j=1:grid.ny],zb2)

@test expected_melt_rate_2d_time1 ≈ basal_melt1
@test expected_melt_rate_2d_time2 ≈ basal_melt2

end


mktempdir() do temporary_directory
    cd(temporary_directory) do
        println("Running climate forcing test in directory:",temporary_directory)
        grid = iceberg_grid()
        climate_forcing_test(grid)
    end
end


