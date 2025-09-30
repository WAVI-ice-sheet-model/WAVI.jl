"""
    JKVstep!(inversion_simulation)

Perform one JKV of the inversion_simulation
"""
function JKVstep!(inversion_simulation)
    
    @unpack model, inversion, JKVstepping_params, output_params, clock = inversion_simulation
   
    update_βeff!(model)
    update_βeff_on_uv_grids!(model)
    update_betas_dirichlet!(model,inversion)
    update_rheological_operators!(model)
    solve_dirichelt_neumann_velocities!(model, inversion, clock)
    update_surf_stress_dirichelt!(inversion)
     

    ###ISSUE is that gh.quad are only on model grid and not also inversion, but I don't want to copy them, or to call model and inversion to all these functions...

    update_velocities_on_h_grid!(model)  
    update_velocities_on_h_grid_dirichlet!(inversion)
    #
    update_shelf_strain_rate!(model)
    update_shelf_strain_rate!(inversion)
    #
    update_av_speed!(model)
    update_av_speed!(inversion)
    #
    update_bed_speed!(model)
    #dirichlet already done above.
    #
    update_surf_speed!(model)
    update_surf_speed!(inversion)
    
    update_surface_velocities_on_uv_grid!(model)
    update_surface_velocities_on_uv_grid!(inversion)

    update_basal_drag_components!(model)
    update_basal_drag_components!(inversion)
    #
    update_dhdt!(model)
    #
    update_shelf_heating!(model)
    update_shelf_heating_dirichlet!(model,inversion)
    #
    update_vert_shear_heating!(model)
    update_vert_shear_heating_dirichlet!(model,inversion)
    #
    update_drag_heating!(model)
    update_drag_heating!(inversion)
   
    update_β_inversion!(model,inversion)
    update_preBfactor_inversion!(model,inversion)
    update_damage!(model)
    update_glen_b!(model)

    inner_update_viscosity!(model)
    update_av_viscosity!(model)
    update_quadrature_falpha!(model)
    update_viscosities_quadratures_dirichlet(model,inversion)
  
    update_JKV!(model,inversion,clock)
    update_JRMS!(model,inversion,clock)
    update_clock_inversion!(inversion_simulation)

    return inversion_simulation
end

"""
    update_clock!(inversion_simulation::AbstractSimulation)

Update the inversion_simulation clock
"""
function update_clock_inversion!(inversion_simulation)
    @unpack clock=inversion_simulation
    clock.n_iter += 1
   # clock.time += JKVstepping_params.dt
    return inversion_simulation
end

"""
    run_inversion_simulation(inversion_simulation)
Perform the inversion_simulation specified by the inversion_simulation
"""
function run_inversion_simulation!(inversion_simulation)
    @unpack model, inversion, JKVstepping_params, output_params = inversion_simulation
    chkpt_tag = "A"



    for   i = (inversion_simulation.clock.n_iter+1): JKVstepping_params.n_iter_total

      if inversion_simulation.clock.n_iter == 0
      update_surface_elevation!(model)
      update_geometry_on_uv_grids!(model)
      update_height_above_floatation!(model)
      update_grounded_fraction_on_huv_grids!(model)
      update_accumulation_rate!(model)
      start_guess_β_inversion!(model,inversion)
      start_guess_η_inversion!(model,inversion)
      update_quadrature_falpha!(model)
      update_av_viscosity!(model)
      update_viscosities_quadratures_dirichlet(model,inversion)
      end

        JKVstep!(inversion_simulation)



        #check if we have hit a temporary checkpoint
        if mod(i,JKVstepping_params.n_iter_chkpt) == 0
            #output a temporary checkpoint
            fname = string("Chkpt",chkpt_tag, ".jld2")
            @save fname inversion_simulation
            chkpt_tag = (chkpt_tag == "A") ? "B" : "A"
            println("making temporary checkpoint at timestep number $(simulation.clock.n_iter)")
        end

        #check if we have hit a permanent checkpoint
        if mod(i,inversion_simulation.JKVstepping_params.n_iter_pchkpt) == 0
            #output a permanent checkpoint
            n_iter_string =  lpad(inversion_simulation.clock.n_iter, 10, "0"); #filename as a string with 10 digits
            fname = string(output_params.output_path, "PChkpt_",n_iter_string, ".jld2")
            @save fname inversion_simulation
            println("making permanent checkpoint at timestep number $(inversion_simulation.clock.n_iter)")
        end

        #check if we have hit an output timestep
        if mod(i,inversion_simulation.output_params.n_iter_out) == 0
            write_output(inversion_simulation)
        end   

    end
    
    #zip the inversion_simulation output (no zipping handled by zip_output)
    zip_output(inversion_simulation)

        
    return inversion_simulation
end