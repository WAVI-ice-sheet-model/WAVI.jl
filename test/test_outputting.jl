using Test, WAVI, MAT, JLD2, NCDatasets, Dates

function output_test(; dt  = 0.5, 
                    end_time = 100., 
                    output_freq = 5., 
                    output_format = "jld2", 
                    zip_format = "none", 
                    prefix = "outfile", 
                    dump_vel = false, 
                    output_path = ".",
                    pchkpt_freq = Inf,
                    chkpt_freq = Inf,
                    niter0 = 0,
                    output_start = false)

    grid = Grid() #default grid with nx = 80, ny = 10
    bed = WAVI.mismip_plus_bed
    solver_params = SolverParams(maxiter_picard = 1)
    model = Model(grid = grid, bed_elevation = bed, solver_params = solver_params)
    timestepping_params = TimesteppingParams(niter0 = niter0, dt = dt, end_time = end_time, chkpt_freq = chkpt_freq, pchkpt_freq = pchkpt_freq)

    outputs = (h = model.fields.gh.h,
                u = model.fields.gh.u,
                v = model.fields.gh.v,
                b = model.fields.gh.b) #output velocities and thickness

    output_freq = output_freq
    output_params = OutputParams(outputs = outputs, 
                        output_freq = output_freq,
                        output_format = output_format,
                        output_path = output_path,
                        zip_format = zip_format,
                        prefix = prefix,
                        dump_vel = dump_vel,
                        output_start = output_start)

    simulation = Simulation(model = model, 
                        timestepping_params = timestepping_params,
                        output_params = output_params)
            
    #perform the simulation
    run_simulation!(simulation)
    return simulation
end

#########################################################
############### Test flags ##############################
test_basic_outputting = true
test_zipping_output   = true
test_output_after_restart = true
test_output_errors = true
##########################################################


@testset "Outputting" begin
    if test_basic_outputting
        @testset "Output files" begin 
            @info "Testing outputting..."
            
            for output_format in ["mat", "jld2"] #check that both mat and jld2 work
                for folder = ["outputs/", "outputs", "./outputs"] #check that both with and without / works
                    isdir(folder) && rm(folder, force = true, recursive = true) #make a new folder
                    mkdir(folder)

                    #run the simulation
                    sim = output_test(output_path = folder, 
                                end_time = 100., 
                                output_freq = 5., 
                                prefix = "testoutfile", 
                                output_format = output_format)

                    foldersim = sim.output_params.output_path 
                    @test foldersim[end] == '/' #test that we do have the / at end of path
                    files = [string(foldersim,f) for f in  readdir(foldersim) if endswith(f, output_format)]
                    println(files)
                    @test length(files) == 20    #check there are the correct number of output files

                    #check they have the correct suffix
                    prefices = [split(str,".")[1] for str in files]
                    suffices = [split(str,".")[end] for str in files]
                    @test all(suffices .== output_format)

                    #check that variables included and correct size
                    fname = string(foldersim,readdir(folder)[1])
                    if output_format == "mat"
                    dict = matread(fname)
                    elseif output_format == "jld2"
                    dict = load(fname)
                    end
                    @test length(dict["t"]) == 1
                    @test size(dict["x"]) == (80,10)
                    @test size(dict["y"]) == (80,10)
                    @test size(dict["b"]) == (80,10)
                    @test size(dict["u"]) == (80,10)
                    @test size(dict["h"]) == (80,10)
                    @test size(dict["v"]) == (80,10)

                    #check it has the correct prefix and number of zeros
                    fname = readdir(folder)[1]
                    @test fname[1:11]  == "testoutfile"

                    #delete the folder
                    rm(folder, force = true, recursive = true)
                end
            end #end loop over output format
        
        end
        @testset "Test output start" begin 
            for output_start in [true, false]
                folder = "outputs/"
                isdir(folder) && rm(folder, force = true, recursive = true) #make a new folder
                mkdir(folder)
                #test that if we only get an output at the start (e.g. if the output freq > end time), then outputted solution matches the expected thing
                sim = output_test(output_path = folder, 
                                end_time = 1., 
                                output_freq = 5., 
                                prefix = "testoutfile", 
                                output_format = "jld2", 
                                output_start = output_start)

                foldersim = sim.output_params.output_path 
                if output_start
                    #check only a single file outputted
                    @test length(readdir(foldersim)) == 1

                    fname = string(foldersim,readdir(folder)[1])
                    dict = load(fname)

                    #check that thickness same as IC and time = 0
                    @test dict["h"] == sim.model.initial_conditions.initial_thickness
                    @test dict["t"] == 0.0


                else
                    #check that no file output
                    @test isempty(readdir(foldersim))
                end

                #delete the folder
                rm(folder, force = true, recursive = true)
            end
        end
    end

    if test_zipping_output
        @testset "Zipping output" begin 
            @info "Testing zipping output..."
            for output_format in ["mat", "jld2"]
                for output_start in [true, false]
                    folder = "outputs/"
                    isdir(folder) && rm(folder, force = true, recursive = true) #make a new folder
                    mkdir(folder)

                    sim = output_test(zip_format = "nc", 
                                        output_path = folder, 
                                        output_format = "mat",
                                        output_start = output_start)
                    
                    @test sim isa Simulation
                    fname = string(folder, sim.output_params.prefix, ".nc")
                    @test isfile(fname) #check the zipped file exists

                    #test variables read from nc file
                    ds = NCDataset(fname)
                    x = ds["x"][:]
                    y = ds["y"][:]
                    h = ds["h"][:,:,:]
                    u = ds["u"][:,:,:]
                    v = ds["v"][:,:,:]
                    b = ds["b"][:,:,:]
                    t = ds["TIME"][:]
                    close(ds)
                    if output_start
                        ds = NCDataset(fname)
                        @test ds["TIME"][:] == 0.:5.:100.
                        close(ds)
                    else
                        ds = NCDataset(fname)
                        @test ds["TIME"][:] == 5.:5.:100.
                        close(ds)
                    end
                        
                    @test size(x) == (80,)
                    @test size(y) == (10,)
                    @test size(b) == (length(x),length(y),length(t))
                    @test size(u) == (length(x),length(y),length(t))
                    @test size(v) == (length(x),length(y),length(t))
                    @test size(h) == (length(x),length(y),length(t))
            
                    #delete the folder
                    rm(folder, force = true, recursive = true)
                end
            end
        end
    end

    if test_output_after_restart
        @testset "Outputting after a restart" begin 
            @info "Testing outputting after a restart"
            
            for output_format in ["mat", "jld2"]

                folder = "outputs/"
                isdir(folder) && rm(folder, force = true, recursive = true) #make a new folder
                mkdir(folder)

                #run simulation first time with no zipping
                output_test(; dt  = 0.5, 
                            end_time = 10., 
                            output_freq = 0.5, 
                            output_format = output_format, 
                            output_path = folder,
                            zip_format = "nc", 
                            dump_vel = true, 
                            pchkpt_freq = 1.)

                #get the time that the first file was outputted (test for https://github.com/RJArthern/WAVI.jl/issues/35)
                first_file_name = joinpath("outputs",string("outfile0000000001.", output_format));
                @test isfile(first_file_name) #check the first output point exists
                dt1 = Dates.unix2datetime(mtime(first_file_name)) #time of last modification

                #run again with different niter0
                sim = output_test(; dt  = 0.5, 
                            end_time = 20., 
                            output_freq = 0.5, 
                            niter0 = 20,
                            output_format = output_format, 
                            output_path = folder,
                            zip_format = "nc", 
                            dump_vel = true, 
                            pchkpt_freq = 1.)

                #check that the time that first file modified has not changed -- i.e. simulation has not touched outputs at earlier times than it should have            
                dt2 = Dates.unix2datetime(mtime(first_file_name)) 
                @test dt1 == dt2
                
                #test variables read from nc file
                fname = string(folder, sim.output_params.prefix, ".nc")
                println(fname)
                ds = NCDataset(fname)
                x = ds["x"][:]
                y = ds["y"][:]
                h = ds["h"][:,:,:]
                u = ds["u"][:,:,:]
                v = ds["v"][:,:,:]
                b = ds["b"][:,:,:]
                t = ds["TIME"][:]
                close(ds)
                ds = NCDataset(fname)
                @test ds["TIME"][:] == 0.5:0.5:20.
                close(ds)
                @test size(x) == (80,)
                @test size(y) == (10,)
                @test size(b) == (length(x),length(y),length(t))
                @test size(u) == (length(x),length(y),length(t))
                @test size(v) == (length(x),length(y),length(t))
                @test size(h) == (length(x),length(y),length(t))

                #check we have the right number of output files 
                foldersim = sim.output_params.output_path 
                @test foldersim[end] == '/' #test that we do have the / at end of path
                files = [string(foldersim,f) for f in  readdir(foldersim) if endswith(f, output_format)]
                println(files)
                @test length(files) == 40    #check there are the correct number of output files
                
                #now repeat for just one timestep
                sim = output_test(; dt  = 0.5, 
                                end_time = 20.5, 
                                output_freq = 0.5, 
                                niter0 = 40,
                                output_format = output_format, 
                                output_path = folder,
                                zip_format = "nc", 
                                dump_vel = true, 
                                pchkpt_freq = 1.)
                
                #test variables read from nc file
                fname = string(folder, sim.output_params.prefix, ".nc")
                println(fname)
                ds = NCDataset(fname)
                x = ds["x"][:]
                y = ds["y"][:]
                h = ds["h"][:,:,:]
                u = ds["u"][:,:,:]
                v = ds["v"][:,:,:]
                b = ds["b"][:,:,:]
                t = ds["TIME"][:]
                close(ds)
                ds = NCDataset(fname)
                @test ds["TIME"][:] == 0.5:0.5:20.5
                close(ds)
                @test size(x) == (80,)
                @test size(y) == (10,)
                @test size(b) == (length(x),length(y),length(t))
                @test size(u) == (length(x),length(y),length(t))
                @test size(v) == (length(x),length(y),length(t))
                @test size(h) == (length(x),length(y),length(t))

                #check we have the right number of files
                foldersim = sim.output_params.output_path 
                @test foldersim[end] == '/' #test that we do have the / at end of path
                files = [string(foldersim,f) for f in  readdir(foldersim) if endswith(f, output_format)]
                println(files)
                @test length(files) == 41    #check there are the correct number of output files

                #check we have dumped the velocity
                @test isfile(string(foldersim, "outfile.nc")) #check the zipped file exists


                # delete everything you just made
                rm(folder, force = true, recursive = true)
                foreach(rm, filter(endswith(".mat"), readdir()))
                foreach(rm, filter(endswith(".jld2"), readdir()))
                foreach(rm, filter(endswith(".bin"), readdir()))
            end
        end
    end

    if test_output_errors
        @testset "Outputting errors" begin 
            @info "Testing outputting errors"
            #check error for weird output format 
            @test_throws ArgumentError output_test(output_format = "incorrect_format")

            #check that non-standard zip_format reverts to none
            sim = output_test(end_time = 1., zip_format = "incorrect_zip_format")
            @test sim.output_params.zip_format == "none"

            @test_throws ArgumentError output_test(output_freq =  -2.)

        end
    end
end

