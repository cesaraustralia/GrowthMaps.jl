
"""
A `RateModel` contains the parameters to calculate the contribution of a
growth or stress factor to the overall population growth rate.

A [`rate`](@ref) method corresponding to the model type must be defined specifying
the formulation, and optionally a [`condition`](@ref) method to e.g. ignore data below
thresholds.

## Returns

An `AbstractFloat` with implicit units of `N*N^-1*T^-1`, where `N` is the 
number of individuals and `T` is time, usually one day. 

It should _not_ return a Unitful.jl `Quantity`.
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
    rate(m::RateModel, x)

Returns a growth rate multiplier for a model given environmental 
variable x, for a single cell. Must be defined for all `RateModel`.

## Arguments

- `m`: the object containing the model parameters
- `x`: the value at a particular location for of an environmental variable
  specified by the models `key()` method, using corresponding to its `key` field.

`rate` methods should use the `@inline` macro to force inlining.
"""
function rate end

"""
    condition(m::RateModel, x)

Subset the data before applying the rate method, given the model `m` and
the environmental variable of value x. 

Returns a `Bool`. `true` will run the `rate` model, `false` will keep the 
current growth rate, which is `1` for a single model.

If not defined for a given `RateModel`, `condition` returns `true` for all cells.

`condition` methods should use the `@inline` macro to force inlining.
"""
function condition end
condition(m, x) = true

abstract type AbstractLowerStress <: StressModel end

@inline rate(m::AbstractLowerStress, x) = (m.threshold - x) * m.mortalityrate
@inline condition(m::AbstractLowerStress, x) = x < m.threshold

"""
    LowerStress(key::Symbol, threshold, mortalityrate)

A [`StressModel`](@ref) where stress occurs below a `threshold` at a
specified `mortalityrate` for some environmental layer `key`.

If units are used, the `Layer` that wraps this models must be in units convertible to the 
units used in `threshold`.
"""
struct LowerStress{T,M} <: AbstractLowerStress
    threshold::T
    mortalityrate::M
end

abstract type AbstractUpperStress <: StressModel end

@inline condition(m::AbstractUpperStress, x) = x > m.threshold
@inline rate(m::AbstractUpperStress, x) = (x - m.threshold) * m.mortalityrate

"""
    UpperStress(key::Symbol, threshold, mortalityrate)

A [`StressModel`](@ref) where stress occurs above a `threshold` at the
specified `mortalityrate`, for some environmental layer `key`.
"""
struct UpperStress{T,M} <: AbstractUpperStress
    threshold::T
    mortalityrate::M
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
struct SchoolfieldIntrinsicGrowth{P,HA,HL,HH,TL,TH,TR} <: GrowthModel
    p::P
    ΔH_A::HA
    ΔH_L::HL
    T_halfL::TL
    ΔH_H::HH
    T_halfH::TH
    T_ref::TR
end

@inline function rate(m::SchoolfieldIntrinsicGrowth, x)
    @fastmath m.p * x/m.T_ref * exp(m.ΔH_A/R * (1/m.T_ref - 1/x)) /
        (1 + exp(m.ΔH_L/R * (1/m.T_halfL - 1/x)) + exp(m.ΔH_H/R * (1/m.T_halfH - 1/x)))
end

# Layers

"""
    Layer{K,U}(model::RateModel)
    Layer(key::Symbol, model::RateModel)
    Layer(key::Symbol, u::Union{Units,Quantity}, model::RateModel)

Layers connect a model to a data source, providing the key to look
up the layer in a GeoData.jl `GeoStack`, and specifying the
scientific units of the layer, if it has units. Using units adds
an extra degree of safety to your calculation, and allows for using
data in different units with the same models.
"""
struct Layer{K,U,M}
    model::M
end
Layer{K,U}(model::RateModel) where {K,U} = Layer{K,U,typeof(model)}(model)
Layer(key::Symbol, model::RateModel) = Layer{key,unit(u"1"),typeof(model)}(model)
Layer(key::Symbol, u::Units, model::RateModel) = Layer{key,u,typeof(model)}(model)
Layer(key::Symbol, q::Quantity, model::RateModel) =
    Layer{key,typeof(q),typeof(model)}(model)

model(l::Layer) = l.model

@inline rate(l::Layer, x) = rate(l.model, x)
@inline condition(l::Layer, x) = condition(l.model, x)

ConstructionBase.constructorof(::Type{<:Layer{K,U}}) where {K,U} = Layer{K,U}

Base.keys(::Layer{K}) where K = K
Base.keys(layers::Tuple{Vararg{<:Layer}}) = map(keys, layers)

Unitful.unit(::Layer{K,U}) where {K,U} = U
