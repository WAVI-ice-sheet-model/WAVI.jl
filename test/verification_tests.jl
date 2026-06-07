using Test, WAVI, LinearAlgebra, CairoMakie, ImageFiltering
@testset  "WAVI tests" begin
    @testset "Iceberg" begin
        @info "Performing tests on a spinning, drifting iceberg."
        include("verification_tests/iceberg_test.jl")
        sim, relerr_h, relerr_u, relerr_v, relerr_theta = iceberg_test(;end_time=1000.)
        @test relerr_h < 1.0e-4
        @test relerr_u < 3.0e-4
        @test relerr_v < 3.0e-4
        @test_broken relerr_theta < 2.0e-4
    end

    @testset "GlaDS" begin
        @info "Performing GlaDS hydrology test on an idealised setup."
        include("verification_tests/GlaDS_test.jl")

        i_shmip_test = 1
        shmip = [7.93e-11, 1.59e-9, 5.79e-9, 2.5e-8, 4.5e-8, 5.79e-7];
        melt_rate = shmip[i_shmip_test]*(3600*24*365.25); # 7.93e-11 is ~2.5 mm/yr, 5.79e-7 is ~18.3 m/yr
        do_visu = false; # create output plot
        sim = GlaDS_test(;dt_days=0.1, end_time_days=5.0, melt_rate=melt_rate, do_visu=do_visu);

        # get reference hydraulic potential and water sheet thickness for chosen shmip melt rate
        include("verification_tests/check-vs-gladsog.jl")
        ref_phi, ref_h = dataOG[shmip[i_shmip_test]];
        x_coords = 1e3 * [5, 15, 25, 35, 45, 55, 65, 75, 85, 95];  # x-coordinates of the reference data points

        # get model indices closest to the reference data coordinates
        grid_x = sim.model.grid.xxh[:,1];  # grid cell center locations of model grid
        calc_h = sim.model.fields.gh.basal_water_thickness;
        y_index = Int(size(calc_h)[2]/2);

        closest_indices = [findmin(abs.(grid_x .- x))[2] for x in x_coords];
        h_err = sum(abs.(calc_h[closest_indices,y_index] .- ref_h[:,1]))/length(x_coords);

        @test h_err < 1e-2 # [1e-3,1e-2,5e-2] for end_time_days=[0.1,5.0,1000]
    end
end

if true
@testset "MISMIP+ verification experiments" begin 
    @testset "MISMIP+ Ice0 verification experiments" begin
        @info "Performing MISMIP+ Ice0 verification experiments: forward and inversion"
        include(joinpath("verification_tests","MISMIP_PLUS_Ice0.jl"))
         include(joinpath("verification_tests","MISMIP_PLUS_inversion.jl"))
        simulation=MISMIP_PLUS_Ice0()
        inversion_simulation=MISMIP_PLUS_inversion(simulation)
        glx=WAVI.get_glx(simulation.model)
        glxtest=glx[[1,div(simulation.model.grid.ny,2),div(simulation.model.grid.ny,2)+1,simulation.model.grid.ny]]
        @test length(glx) == simulation.model.grid.ny #check that the grounding line covers the whole domain in the y-direction
        @test (glxtest[4]-glxtest[1])/(glxtest[4]+glxtest[1]) < 1e-4
        @test (glxtest[2]-glxtest[3])/(glxtest[2]+glxtest[3]) < 1e-4
        @testset "Tight Tolerance" begin
            # For quick testing at low resolutions (e.g. 8 km) these may be broken
            @test_broken 480000<glxtest[1]<540000
            @test_broken 480000<glxtest[4]<540000
            @test_broken 430000<glxtest[2]<460000
            @test_broken 430000<glxtest[3]<460000
        end
        @testset "Loose Tolerance" begin
            @test 480000<glxtest[1]<550000
            @test 480000<glxtest[4]<550000
            @test 430000<glxtest[2]<490000
            @test 430000<glxtest[3]<490000
        end
        @test all((inversion_simulation.model.fields.gh.preBfactor) .< 1.1)
        @test all((inversion_simulation.model.fields.gh.preBfactor) .> 0.9)

        #check the melt rate for ice_1r is doing something sensible
        function m1(h, b)
            draft = -(918.0 / 1028.0) * h
            cavity_thickness = draft .- b
            cavity_thickness = max.(cavity_thickness, 0)
            m =  0.2*tanh.(cavity_thickness./75).*max.((-100 .- draft), 0)
            return m
        end
        melt = m1.(simulation.model.fields.gh.h, simulation.model.fields.gh.h)
        @test all(melt .>= 0)
        @test maximum(melt) < 100

    end
end
end