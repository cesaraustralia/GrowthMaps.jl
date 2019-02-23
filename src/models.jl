
abstract type AbstractRateModel end

abstract type AbstractGrowthModel <: AbstractRateModel end
abstract type AbstractStressModel <: AbstractRateModel end


scale(model::AbstractRateModel) = model.scale

"""
Let the framework know which data to supply, in which units
"""
function withdata end

"""
For subseting the matrix before applying the main rule
Returns a boolean
"""
function subset end

"""
For calculative the cumulative rate for a single cell
"""
function rule end

datakey(x) = withdata(x)[1]
dataunit(x) = withdata(x)[2]


struct DryDays{L,S} <: AbstractStressModel
    lowersm::L
    scale::S
end

@inline withdata(::DryDays) = :sm_surface, 1
@inline subset(m::DryDays, sm) = sm < m.lowersm
@inline rule(m::DryDays, x) = abs(x - m.lowersm)


struct ColdDays{L,S} <: AbstractStressModel
    lowerct::L
    scale::S
end

@inline withdata(::ColdDays) = :surface_temp, K
@inline subset(m::ColdDays, temp) = temp < m.lowerct
@inline rule(m::ColdDays, x) = abs(x - m.lowerct)

 

struct HotDays{U,S} <: AbstractStressModel
    upperct::U
    scale::S
end

@inline withdata(::HotDays) = :surface_temp, K
@inline subset(m::HotDays, temp) = temp > m.upperct
@inline rule(m::HotDays, x) = x - m.upperct



struct Wilting{W,S} <: AbstractStressModel
    lowerwilt::W
    scale::S
end

@inline withdata(::Wilting) = :land_fraction_wilting, 1
@inline subset(m::Wilting, wilting) = wilting > m.lowerwilt
@inline rule(m::Wilting, x) = x - m.lowerwilt



struct IntrinsicPopGrowth{P,K,S} <: AbstractGrowthModel
    p25::P
    H_A::K
    H_L::K
    T_0_5L::K
    H_H::K
    T_0_5H::K
    scale::S
end

const R = 1.987

@inline withdata(::IntrinsicPopGrowth) = :surface_temp, K
@inline subset(m::IntrinsicPopGrowth, temp) = true
@inline rule(m::IntrinsicPopGrowth, x) = m.p25 * (x * (m.H_A/R * (1/K(25.0°C) - 1/x))) /
    (1 + exp(m.H_L/R * (1/m.T_0_5L - 1/x)) + exp(m.H_H/R * (1/m.T_0_5H - 1/x)))
