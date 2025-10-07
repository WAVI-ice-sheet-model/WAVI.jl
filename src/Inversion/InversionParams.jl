struct InversionParams{T <: Real}
                    reltol :: T
                    abstol :: T
                    maxiter :: Int
                    gmres_restart :: Int
                    βgrounded_start :: T
                    βfloating_start :: T
                    ηstart_guess :: T
                    βpower :: T
                    Bpower_shelf :: T
                    Bpower_grounded :: T
                    inner_tol ::T 
                    inner_maxiters:: Int
                    cg:: Bool
                    gmres:: Bool
end

"""
InversionParams(; <kwargs>)

Construct a WAVI.jl parameters object for holding physical parameters.

Keyword arguments
=================
- `reltol`: relative tolerance 
- `abstol`: asbolute tolerance
- `maxiter`: maximum number of iterations in matrix solves 
- ` gmres_restart`: number of iterations before restarting in gmres
- `βgrounded_start`: initial guess of β over grounded ice (scalar)
- `βfloating_start`: initial guess of β over floating ice (scalar)
- `ηstart_guess`: initial guess of η (scalar)
- `βpower`: p in equation (36) in Arthern et al, 2015
- `Bpower_shelf`: p in equation (36) in Arthern et al, 2015
- `Bpower_grounded`: p in equation (38) in Arthern et al, 2015
- `inner_tol`: tolerance for the inner problem
- `inner_maxiters`: maximum number of iterations for the inner problem
- `cg`: use congugate gradient
- `gmres`: use gmres 
"""

function InversionParams(; reltol::T = 1e-6, 
                    abstol::T = 1e-6, 
                    maxiter::Int = 1000,
                    gmres_restart::Int=50,
                    βgrounded_start::T=1.e4,
                    βfloating_start::T=1.e-4,
                    ηstart_guess::T = 1.0e7,
                    βpower::T = 0.1,
                    Bpower_shelf::T = 0.1,
                    Bpower_grounded::T = 0.01,
                    inner_tol::T = 1.e-4,
                    inner_maxiters::Int=1000,
                    cg::Bool=false,
                    gmres::Bool=true) where {T <: Real}
                      
  return InversionParams{T}(reltol,
                  abstol, 
                  maxiter,
                  gmres_restart,
                  βgrounded_start,
                  βfloating_start,
                  ηstart_guess,
                  βpower,
                  Bpower_shelf,
                  Bpower_grounded,
                  inner_tol, 
                  inner_maxiters,
                  cg,
                  gmres
                  )
end
