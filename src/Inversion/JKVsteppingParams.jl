
struct JKVsteppingParams{N <: Integer, T <: Real}
            niter0 :: N       #starting iteration number
      n_iter_chkpt :: T       #number of iterations per temporary checkpoint
     n_iter_pchkpt :: T       #number of iterations per permanent checkpoint 
        n_iter_out :: T       #number of iterations per output
       chkpt_path  :: String  #path to location of permanent and tempoary checkpoint output
      n_iter_total :: N       #total number of JKV iterations counting from zero
   end


"""
Construct a WAVI.jl JKVsteppingParams object.
JKVsteppingParams stores information relating to JKV stepping in inversion.

Keyword arguments
=================
- niter0:                   Iteration number of the first JKVstep. niter0 = 0 corresponds to a new simulation, while niter0 > 0 (positive integer) corresponds to a pickup.
- n_iter_chkpt:             Frequency of outputting temporary checkpoints
- n_iter_pchkpt:            Frequency of outputting permenant checkpoints
- n_iter_out:               Number of iterations per output
- chkpt_path:               Path to location checkpoint output
- n_iter_total:             Total number of JKV iterations counting from zero
"""

function JKVsteppingParams(;
                        niter0 = 0,
                        n_iter_chkpt = Inf,
                        n_iter_pchkpt = Inf,
                        n_iter_out = Inf,
                        chkpt_path = "./",
                        n_iter_total = 0)

    #check the output path ends in '/' and exists
    endswith(chkpt_path, "/") || (chkpt_path = string(chkpt_path, "/"))
    if ~isdir(chkpt_path)
        @warn string("Did not find output path ", chkpt_path, ". Any outputs will go to the working directory", pwd())
        chkpt_path = "./"
    end

    return JKVsteppingParams(niter0,n_iter_chkpt, n_iter_pchkpt,n_iter_out, 
                            chkpt_path,n_iter_total)
end
