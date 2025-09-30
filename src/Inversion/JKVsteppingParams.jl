
struct JKVsteppingParams{N <: Integer}
            niter0 :: N       #starting iteration number
      n_iter_chkpt :: N       #number of iterations per temporary checkpoint
     n_iter_pchkpt :: N       #number of iterations per permanent checkpoint 
       chkpt_path  :: String  #path to location of permanent and tempoary checkpoint output
      n_iter_total :: N       #total number of JKV iterations counting from zero
        n_iter_out :: N       #number of iterations per output
   end


"""
Construct a WAVI.jl JKVsteppingParams object.
JKVsteppingParams stores information relating to JKV stepping in inversion.

Keyword arguments
=================
- niter0:                   Iteration number of the first JKVstep. niter0 = 0 corresponds to a new simulation, while niter0 > 0 (positive integer) corresponds to a pickup.
- n_iter_chkpt:             Frequency of outputting temporary checkpoints
- n_iter_pchkpt:            Frequency of outputting permenant checkpoints
- chkpt_path:               Path to location checkpoint output
- n_iter_total:             Total number of JKV iterations counting from zero
- n_iter_out:               Number of iterations per output
"""

function JKVsteppingParams(;
                        niter0 = 0,
                        n_iter_chkpt = Inf,
                        n_iter_pchkpt = Inf,
                        chkpt_path = "./",
                        n_iter_total = 5, 
                        n_iter_out = 5)

    #check the output path ends in '/' and exists
    endswith(chkpt_path, "/") || (chkpt_path = string(chkpt_path, "/"))
    if ~isdir(chkpt_path)
        @warn string("Did not find output path ", chkpt_path, ". Any outputs will go to the working directory", pwd())
        chkpt_path = "./"
    end

    return JKVsteppingParams(niter0,n_iter_chkpt, n_iter_pchkpt, 
                            chkpt_path,n_iter_total, n_iter_out)
end
