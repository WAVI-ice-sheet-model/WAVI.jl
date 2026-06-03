export ISMIP7

struct ISMIP7{S <: String} <: AbstractClimateForcing
   gcm::S
   scenario::S
end

ISMIP7(;gcm = "CESM2-WACCM",scenario = "585") = ISMIP7(gcm,scenario)