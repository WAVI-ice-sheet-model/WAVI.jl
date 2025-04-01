struct InversionParams{T <: Real, A, W, G}
                    gmres_reltol :: T
                    gmres_abstol :: T
                    gmres_maxiter :: Int
                    gmres_restart :: Int
                    accumulation_rate_holder :: A
                    glen_a_ref_holder :: G
                    weertman_c_holder :: W
                    βgrounded_start :: T
                    βfloating_start :: T
                    ηstart_guess:: T
                    βpower :: T
                    Bpower_shelf :: T
                    Bpower_grounded :: T
                    inner_tol ::T 
                    inner_maxiters:: T
end

"""
InversionParams(; <kwargs>)

Construct a WAVI.jl parameters object for holding physical parameters.

Keyword arguments
=================
- 
"""
function InversionParams(; gmres_reltol = 1e-6, 
                    gmres_abstol = 1e-6, 
                    gmres_maxiter = 1000,
                    gmres_restart=100,
                    accumulation_rate_holder=0.0,
                    glen_a_ref_holder=0.0,
                    weertman_c_holder=0.0,
                    βgrounded_start=1.e4,
                    βfloating_start=1.e-4,
                    ηstart_guess = 1.0e7,
                    βpower = 0.1,
                    Bpower_shelf = 0.1,
                    Bpower_grounded = 0.01,
                    inner_tol = 1.e-4,
                    inner_maxiters=1000 )
                      
  return InversionParams(gmres_reltol,
                  gmres_abstol, 
                  gmres_maxiter,
                  gmres_restart,
                  accumulation_rate_holder,
                  glen_a_ref_holder,
                  weertman_c_holder,
                  βgrounded_start,
                  βfloating_start,
                  ηstart_guess,
                  βpower,
                  Bpower_shelf,
                  Bpower_grounded,
                  inner_tol, 
                  inner_maxiters
                  )
end
