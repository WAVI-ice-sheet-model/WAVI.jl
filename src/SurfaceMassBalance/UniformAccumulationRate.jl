export UniformAccumulationRate

struct UniformAccumulationRate{T <: Real} <: AbstractSurfaceMassBalance
    a :: T #uniform accumulation rate applied everywhere
end

UniformAccumulationRate(; a = 0.0) = UniformAccumulationRate(a) 

"""
    update_accumulation_rate!(accumulation_rate::UniformAccumulationRate, model::AbstractModel, clock) 

Update the accumulation rate for the UniformAccumulationRate type

"""
function update_accumulation_rate!(surface_mass_balance::UniformAccumulationRate, model::AbstractModel, clock::Clock) 
    @unpack accumulation = model.fields.gh
    accumulation .= surface_mass_balance.a
    return nothing
end