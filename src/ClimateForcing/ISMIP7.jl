export ISMIP7_ANOMALY,ISMIP7_OCX,ISMIP7_CONTROL


struct ISMIP7_ANOMALY{S <: String} <: AbstractClimateForcing
   gcm::S
   scenario::S
end

ISMIP7_ANOMALY(;gcm = "CESM2-WACCM",scenario = "585") = ISMIP7_ANOMALY(gcm,scenario)

struct ISMIP7_CONTROL{S <: String} <: AbstractClimateForcing
   gcm::S
end

ISMIP7_CONTROL(;gcm = "CESM2-WACCM") = ISMIP7_CONTROL(gcm)

struct ISMIP7_OCX <: AbstractClimateForcing end