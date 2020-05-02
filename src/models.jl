"""
Passed to an optimiser to facilitate fitting any `RateModel`,
without boilerplate or methods rewrites.
"""
mutable struct ModelWrapper{M}
    model::M
end
ModelWrapper(m::Tuple) = ModelWrapper{typeof(m)}(m)
ModelWrapper(m...) = ModelWrapper{typeof(m)}(m)

struct Layer{K,M}
    model::M
end
Layer{K}(x) where K = Layer{K,typeof(x)}(x)
Layer(key::Symbol, x) = Layer{key,typeof(x)}(x)

ConstructionBase.constructorof(::Type{<:Layer{K}}) where K = Layer{K}

Base.keys(model::Layer{K}) where K = K
Base.keys(models::Tuple{Vararg{<:Layer}}) = map(keys, models)

"""
Update the model with passed in params, and run [`rate`] over the independent variables.
"""
(f::ModelWrapper)(xs, params) =
    GrowthMaps.conditionalrate.(Ref(reconstruct(f.model, params, Real)), xs)

"""
A `RateModel` conatains the parameters to calculate the contribution of a
growth or stress factor to the overall population growth rate.

A [`rate`](@ref) method corresponding to the model type must be defined specifying
the formulation, and optionally a [`condition`](@ref) method to e.g. ignore data below
thresholds.
"""
abstract type RateModel end

"""
The intrinsic rate of population growth is the exponential growth rate of a
population when growth is not limited by density dependent factors.

More formally, if the change in population size ``N`` with time ``t`` is
expressed as ``\\frac{dN}{dt} = rN``, then ``r`` is the intrinsic population
growth rate (individuals per day per individual) which can be decomposed into
per capita reproduction and mortality rate. The intrinsic growth rate parameter
``r`` depends strongly on temperature, with population growth inhibited at low
and high temperatures (Haghani et al. 2006).

This can be described using a variety of non-linear functions,
see [`SchoolfieldIntrinsicGrowth`](@ref) for an implementation.
"""
abstract type GrowthModel <: RateModel end

"""
Extreme stressor mortality can be assumed to occur once an environmental
variable ``s`` exceeds some threshold (e.g. critical thermal maximum), beyond
which the mortality rate scales approximately linearly with the depth of the
stressor (Enriquez and Colinet 2017).

Stressor induced mortality is incorporated in growth rates by quantifying the
`threshold` parameter ``s_c`` beyond which stress associated mortality
commences, and the `mortalityrate` parameter ``m_s`` which reflects the
per-capita mortality per stress unit, per time (e.g. degrees beyond the stress
threshold per day).

When stressors are combined in a model, the mortality rate for each stressor
``s`` is incorporated as ``\\frac{dN}{dt}=(r_p-r_n)N`` where ``r_n = \\sum{s}f(s,
c_s)m_s`` and ``f(s, c_s)`` is a function that provides the positive units by
which s exceeds ``s_c``.

This relies on the simplifying assumption that different stressors contribute
additively to mortality rate, but conveniently allows the mortality response of
stressors to be partitioned and calculated separately to the intrinsic growth
rate, which aids parameterisation and modularity. More complicated stress
responses are more likely to capture reality more completely, but the
availability of data on detailed mortality responses to climatic stressors is
often lacking.

Upper and lower stresses may be in relation to any environmental variable,
specied with the `key` parameter, that will use the data at that key in `GeoStack`
for the current timestep.
"""
abstract type StressModel <: RateModel end

"""
    rate(m, x)
Calculates the rate modifyer for a single cell. Must be defined for all models.

`m` is the object containing the model parameters, and
`x` is the value at a particular location for of an environmental variable
specified by the models `key()` method, using corresponding to its `key` field.
"""
function rate end
rate(l::Layer, x) = rate(l.model, x)

"""
    condition(m, x)

Subseting the data before applying the rate method, given the model `m` and
the environmental variable of value x. Returns a boolean. If not defined, defaults
to `true` for all values.
"""
function condition end
condition(l::Layer, x) = condition(l.model, x)
condition(m, x) = true

abstract type AbstractLowerStress <: StressModel end

# TODO set units in the model, this is a temporary hack
@inline condition(m::AbstractLowerStress, x) = condition(m, x * oneunit(m.threshold))
@inline condition(m::AbstractLowerStress, x::Quantity) = x < m.threshold
@inline rate(m::AbstractLowerStress, x::Quantity) =
    (m.threshold - x) * m.mortalityrate
@inline rate(m::AbstractLowerStress, x) = rate(m, x * oneunit(m.threshold))

"""
    LowerStress(key::Symbol, threshold, mortalityrate)

A [`StressModel`](@ref) where stress occurs below a `threshold` at a
specified `mortalityrate` for some environmental layer `key`.

Independent variables must be in the same units as
"""
@bounds struct LowerStress{T,M} <: AbstractLowerStress
    threshold::T     | (250, 350)
    mortalityrate::M | (-0.5, 0.0)
end

@bounds struct ColdStress{T,M} <: AbstractLowerStress
    threshold::T     | (250, 350)
    mortalityrate::M | (-0.5, 0.0)
end

@bounds struct MoistureStress{T,M} <: AbstractLowerStress
    threshold::T     | (0.0, 1.0)
    mortalityrate::M | (-0.5, 0.0)
end

abstract type AbstractUpperStress <: StressModel end

@inline condition(m::AbstractUpperStress, x::Quantity) = x > m.threshold
@inline condition(m::AbstractUpperStress, x) = condition(m, x * oneunit(m.threshold))
@inline rate(m::AbstractUpperStress, x::Quantity) =
    (x - m.threshold) * m.mortalityrate
@inline rate(m::AbstractUpperStress, x) = rate(m, x * oneunit(m.threshold))

"""
    UpperStress(key::Symbol, threshold, mortalityrate)

A [`StressModel`](@ref) where stress occurs above a `threshold` at the
specified `mortalityrate`, for some environmental layer `key`.
"""
@bounds struct UpperStress{T,M} <: AbstractUpperStress
    threshold::T     | (250, 350)
    mortalityrate::M | (-0.5, 0.0)
end

@bounds struct HeatStress{T,M} <: AbstractUpperStress
    threshold::T     | (250, 350)
    mortalityrate::M | (-0.5, 0.0)
end

@bounds struct WiltStress{T,M} <: AbstractUpperStress
    threshold::T     | (0.0, 1.0)
    mortalityrate::M | (-0.5, 0.0)
end

"""
SchoolfieldIntrinsicGrowth(p, ΔH_A, ΔH_L, ΔH_H, Thalf_L, Thalf_H, T_ref, R)

A [`GrowthModel`](@ref) where the temperature response of positive growth rate is
modelled following Schoolfield et al. 1981, _"Non-linear regression of biological
temperature-dependent rate models base on absolute reaction-rate theory"_.

The value of the specified data layer _must_ be in Kelvin.

## Arguments
- `p::P`: growth rate at reference temperature `T_ref`
- `ΔH_A`: enthalpy of activation of the reaction that is catalyzed by the enzyme, in cal/mol or J/mol.
- `ΔH_L`: change in enthalpy associated with low temperature inactivation of the enzyme.
- `ΔH_H`: change in enthalpy associated with high temperature inactivation of the enzyme.
- `T_halfL`: temperature (Unitful K) at which the enzyme is 1/2 active and 1/2 low temperature inactive.
- `T_halfH`: temperature (Unitful K) at which the enzyme is 1/2 active and 1/2 high temperature inactive.
- `T_ref`: Reference temperature in kelvin.

This can be parameterised from empirical data, (see Zhang et al. 2000; Haghani
et al. 2006; Chien and Chang 2007) and non-linear least squares regression.
"""
@bounds @flattenable struct SchoolfieldIntrinsicGrowth{P,HA,HL,HH,TL,TH,TR} <: GrowthModel
    p::P         | true  | (0.0, 1.0)
    ΔH_A::HA     | true  | (0.0u"cal/mol", 1e6u"cal/mol")
    ΔH_L::HL     | true  | (-1e6u"cal/mol", 0.0u"cal/mol")
    T_halfL::TL  | true  | (150.0K, 400K)
    ΔH_H::HH     | true  | (0.0u"cal/mol", 1e6u"cal/mol")
    T_halfH::TH  | true  | (150.0K, 400K)
    T_ref::TR    | false | _
end

@inline rate(m::SchoolfieldIntrinsicGrowth, x::Real) =
    rate(m::SchoolfieldIntrinsicGrowth, x * u"K")
@inline rate(m::SchoolfieldIntrinsicGrowth, x::Quantity) = begin
    @fastmath m.p * x/m.T_ref * exp(m.ΔH_A/R * (1/m.T_ref - 1/x)) /
        (1 + exp(m.ΔH_L/R * (1/m.T_halfL - 1/x)) + exp(m.ΔH_H/R * (1/m.T_halfH - 1/x)))
end
