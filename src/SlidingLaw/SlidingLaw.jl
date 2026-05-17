module SlidingLaw

using Parameters

export update_β_using_sliding_law!

using WAVI: AbstractSlidingLaw, AbstractModel


#add each of the individual sliding laws
include("./WeertmanSlidingLaw.jl")
include("./coulomb.jl")
include("./budd.jl")
include("./tsai.jl")
include("./tsaiBudd.jl")
include("./schoof.jl")
include("./zoetIverson.jl")



end