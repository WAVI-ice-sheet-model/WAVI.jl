module SlidingLaw

using Parameters

export update_β_using_sliding_law!

using WAVI: AbstractSlidingLaw, AbstractModel
using WAVI.Grids: Grid


#add each of the individual sliding laws
include("./WeertmanSlidingLaw.jl")
include("./coulomb.jl")
include("./budd.jl")
include("./tsai.jl")
include("./tsaiBudd.jl")
include("./schoof.jl")
include("./zoetIverson.jl")


"""
    update_drag_coefficient!(model::AbstractModel)

Update drag coefficient used in the sliding law to account for migration of grounding line.
This is currently only used for the Weertman sliding law and Weertman part of the Tsai sliding law, 
as they are the only sliding laws that do not directly depend on effective pressure, which already
acounts for migration of grounding line.
"""
function update_drag_coefficient!(model::AbstractModel)
    @unpack gh=model.fields
    @unpack sliding_law=model
    gh.drag_coefficient .= sliding_law.drag_coefficient .* gh.grounded_fraction
    return model
end

function reconstruct_on_grid(sliding_law::SL, grid::Grid) where {SL <: AbstractSlidingLaw}
    return sliding_law
end

function reconstruct_on_subdomain(sliding_law::SL, grid::Grid,subdomain::NTuple{4,<: Integer}) where {SL <: AbstractSlidingLaw}
    return sliding_law
end

end