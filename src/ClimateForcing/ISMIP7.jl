export ISMIP7

struct ISMIP7{S <: String} <: AbstractClimateForcing
   gcm::S
   scenario::S
end

ISMIP7(;gcm = "cesm2waccm",scenario = "ssp585") = ISMIP7(gcm,scenario)