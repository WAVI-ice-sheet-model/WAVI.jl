
struct JKVsteppingParams{N <: Integer}
            niter0 :: N      #starting iteration number
      n_iter_chkpt :: N       #number of iterations per temporary checkpoint
     n_iter_pchkpt :: N      #number of iterations per permanent checkpoint 
       chkpt_path  :: String #path to location of permanent and tempoary checkpoint output
      n_iter_total :: N     #total number of JKVsteps counting from zero
        n_iter_out :: N    #number of iterations per output
   end

function JKVsteppingParams(;
                        niter0 = 0,
                        n_iter_chkpt = 0,
                        n_iter_pchkpt = 0,
                        chkpt_path = "./",
                        n_iter_total = 5, 
                        n_iter_out = 5)

   #=  #check compatibility of n_iter_total and end_time, and compute them 
    end_time, n_iter_total = compute_iterations_and_end_time(end_time, n_iter_total,dt)

    #compute number of timesteps checkpoint number of timesteps
    chkpt_freq == Inf ? n_iter_chkpt = Inf : n_iter_chkpt  = round(Int, chkpt_freq/dt)
    pchkpt_freq == Inf ? n_iter_pchkpt = Inf : n_iter_pchkpt = round(Int, pchkpt_freq/dt)
     =#

    #check the output path ends in '/' and exists
    endswith(chkpt_path, "/") || (chkpt_path = string(chkpt_path, "/"))
    if ~isdir(chkpt_path)
        @warn string("Did not find output path ", chkpt_path, ". Any outputs will go to the working directory", pwd())
        chkpt_path = "./"
    end

    return JKVsteppingParams(niter0,n_iter_chkpt, n_iter_pchkpt, 
                            chkpt_path,n_iter_total, n_iter_out)
end
