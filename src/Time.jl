module Time

import Base: show

export Clock

#mutable clock structure to store time info
mutable struct Clock{T <: Real, N <: Integer}
    n_iter::N
    time::T
end

#clock constructor
function Clock(;
                n_iter = 0,
                time = 0)
    return Clock(n_iter, time)
end

"""
    compute_iterations_and_end_time(end_time, n_iter_total, dt)

Check input end time and total iterations, and compute those not passed.
"""
function compute_iterations_and_end_time(end_time, n_iter_total, dt)
    (isa(dt, Real) && (dt > 0))  || throw(ArgumentError("timestep dt must be a positive number"))

    if (~(n_iter_total === nothing) && ~(end_time === nothing)) #if both passed, throw error if incompatible
        (end_time ≈ (n_iter_total * dt)) ||  throw(ArgumentError("You have specified both end time (end_time) and total iterations (n_iter_total), but their values are incompatible: end time mustequal n_iter_total * dt"))
    elseif ((n_iter_total === nothing) && ~(end_time === nothing)) #if only end time passed, n_iter_total is the nearest integer
        end_time == Inf ? n_iter_total = Inf : n_iter_total  = round(Int, end_time/dt)
    elseif (~(n_iter_total === nothing) && (end_time === nothing)) #if only number of iterations
        n_iter_total == Inf ? end_time = Inf : end_time  = dt*n_iter_total
    elseif (n_iter_total === nothing) && (end_time === nothing) #if neither is passed
        throw(ArgumentError("You must pass at least one of end_time or n_iter_total"))
    end 

    #check both are positive numbers, now that both end_time and n_iter_total assigned
    (isa(n_iter_total, Integer) && (n_iter_total > 0)) || throw(ArgumentError("n_iter_total must be a positive integer"))
    (isa(end_time, Real) && (end_time > 0))|| throw(ArgumentError("end_time must be a positive number"))
    return end_time, n_iter_total
end

function Base.show(io::IO, clock::Clock)
    return print(io, "$(clock.time) iteration $(clock.n_iter)")
end

end
