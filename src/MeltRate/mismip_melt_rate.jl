struct MISMIPMeltRateOne{T <: Real} <: AbstractMeltRate 
    α  :: T
    ρi :: T
    ρw :: T
end

"""
    function MISMIPMeltRateOne(; <kwargs>)

Construct a melt rate to specify the melt rate according to the MISMIP_1r experiment

Keyword arguments
=================
- α : Calibrated melt rate prefactor (defaults to unity, value used in MISMIP+)
- ρi : Ice density
- ρw : Water density
"""
MISMIPMeltRateOne(; α = 1.0, ρi = 918.0, ρw = 1028.0) = MISMIPMeltRateOne(α,ρi, ρw)

"""
    update_shelf_melt_rate(shelf_melt_rate::MISMIPMeltRateOne, fields, grid, clock) 

Update the melt rate under ice shelves for the MISMIPMeltRateOne type, used in MISMIP Ice 1x experiments
"""
function update_shelf_melt_rate!(shelf_melt_rate::MISMIPMeltRateOne, fields, grid, clock) 
    draft = -(shelf_melt_rate.ρi / shelf_melt_rate.ρw) .* fields.gh.h
    cavity_thickness = draft .- fields.gh.b
    cavity_thickness = max.(cavity_thickness, 0)
    m =  shelf_melt_rate.α .* 0.2*tanh.(cavity_thickness./75).*max.((-100 .- draft), 0)
    fields.gh.shelf_basal_melt[:] .= m[:]
end

