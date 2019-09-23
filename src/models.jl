
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
@inline condition(m::LowerStress, x) = x * oneunit(m.threshold) < m.threshold
@inline rate(m::LowerStress, x) = (m.threshold - x * oneunit(m.threshold)) * m.mortalityrate

struct UpperStress{T,M} <: StressModel
    key::Symbol
    threshold::T
    mortalityrate::M
end
@inline condition(m::UpperStress, x) = x * oneunit(m.threshold) > m.threshold
@inline rate(m::UpperStress, x) = (x * oneunit(m.threshold) - m.threshold) * m.mortalityrate

"""
See Schoolfield 1981 "Non-linear regression of biological temperature-dependent rate
models base on absolute reaction-rate theory"
"""
struct IntrinsicGrowth{P,CM,K,R} <: GrowthModel
    key::Symbol
    p::P
    ΔH_A::CM
    ΔH_L::CM
    ΔH_H::CM
    T_halfL::K
    T_halfH::K
    T_ref::K
    R::R
end

rate(m::IntrinsicGrowth, x) = begin
    u = unit(m.T_ref)
    x *= u 
    m.p * x/m.T_ref * exp(m.ΔH_A/m.R * (1/m.T_ref - 1/x)) /
        (1 + exp(m.ΔH_L/m.R * (1/m.T_halfL - 1/x)) + exp(m.ΔH_H/m.R * (1/m.T_halfH - 1/x)))
end

