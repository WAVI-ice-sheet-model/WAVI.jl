export ConstantDamage

struct ConstantDamage <: AbstractFracture end

function update_damage!(fracture::ConstantDamage,model::AbstractModel{T,N};kwargs...) where {T,N}
  return model
end

function update_strain_history!(fracture::ConstantDamage,model::AbstractModel{T,N};kwargs...) where {T,N}
  return model
end