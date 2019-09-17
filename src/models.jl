
abstract type RateModel end

abstract type GrowthModel <: RateModel end
abstract type StressModel <: RateModel end

Base.keys(model::RateModel) = model.key
Base.keys(models::Tuple{Vararg{<:RateModel}}) = tuple(union(keys.(models))...)

"""
For calculative the cumulative rate for a single cell.
"""
function rate end

"""
For subseting the data before applying the rate method.
Returns a boolean
"""
function condition end
condition(m, x) = true


struct LowerStress{T,M} <: StressModel
    key::Symbol
    threshold::T
    mortalityrate::M
end
@inline condition(m::LowerStress, x) = x < m.lowerct
@inline rate(m::LowerStress, x) = (m.lowerct - x) * m.lowerctm 
 
struct UpperStress{T,M} <: StressModel
    key::Symbol
    threshold::T
    mortalityrate::M
end
@inline condition(m::UpperStress, x) = x > m.upperct
@inline rate(m::UpperStress, x) = (x - m.upperct) * m.upperctm

"""
See Schoolfield 1981 "Non-linear regression of biological temperature-dependent rate 
models base on absolute reaction-rate theory"
"""
struct IntrinsicGrowth{S,P,K,TR} <: GrowthModel
    key::Symbol # key::Sym
    scale::S    # scale::S
    p25::P      # p25::P
    H_A::K      # H_A::K
    H_L::K      # H_L::K
    T_0_5L::K   # T_0_5L::
    H_H::K      # H_H::K
    T_0_5H::K   # T_0_5H::
    tref::TR    # tref::TR
end

const R = 1.987
@inline rate(m::IntrinsicPopGrowth, x) = m.p25 * (x * (m.H_A/R * (1/m.tref - 1/x))) /
    (1 + exp(m.H_L/R * (1/m.T_0_5L - 1/x)) + exp(m.H_H/R * (1/m.T_0_5H - 1/x))) / m.scale

