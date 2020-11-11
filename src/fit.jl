
"""
    fit(model, obs::AbstractArray)

Fit a model to data with least squares regression, using `curve_fit` from
LsqFit.jl. The passed in model should be initialised with sensible defaults,
these will be used as the initial parameters for the optimization.

Any (nested) `Real` fields on the struct are flattened to a parameter vector using
[Flatten.jl](http://github.com/rafaqz/Flatten.jl). Fields can be marked to ignore
using the `@flattenable` macro from [FieldMetadata.jl](http://github.com/rafaqz/FieldMetadata.jl).

## Arguments

- `model`: Any constructed [`RateModel`](@ref) or a `Tuple` of `RateModel`.
- `obs`: A `Vector` of `(val, rate)` tuples where `val` is the value of the
  x-axis variable (such as temperature), and `rate` is the growth rate observed.

## Returns

An updated `Model` containing the fitted parameter values.
"""
function fit!(model::Model, obs::AbstractArray)
    fit = curve_fit(first.(obs), last.(obs), collect(model[:val])) do xs, vals
        model[:val] = vals
        GrowthMaps.conditionalrate.(Ref(stripparams(model)), xs)
    end
    model[:val] = fit.param
    return model
end

"""
    manualfit!(model::Model, data::NamedTuple; obs=[], kw...) =

Returns the fitted model.

# Arguments

- `obs`: A `Vector` of `(val, rate)` tuples where `val` is the value of the
  x-axis variable (such as temperature), and `rate` is the growth rate observed.
- `data`: A `NamedTuple` of `AbstractVector`. The `NamedTuple` key/s must match
  the key required by the `Layer`/s.
- `obs`: Optional observations to scatter-plot over the curve. A `Vector` of `(val, rate)`
  tuples where `val` is the value of the x-axis variable (such as temperature), and `rate`
  is the growth rate observed.
- `kwargs`: passed to `plot`

# Example

```julia
p = 3e-01
ΔH_A = 3e4cal/mol
ΔH_L = -1e5cal/mol
ΔH_H = 3e5cal/mol
Thalf_L = 2e2K
Thalf_H = 3e2K
T_ref = K(25.0°C)
growthmodel = SchoolfieldIntrinsicGrowth(p, ΔH_A, ΔH_L, Thalf_L, ΔH_H, Thalf_H, T_ref)
model = Model(Layer(:surface_temp, K, growthmodel))
obs = []

tempdata=(surface_temp=(270.0:0.1:310.0)K,)
manualfit!(model, tempdata; obs=obs)

To use the interface in a desktop app, use Blink.jl:

```julia; eval=false
using Blink
w = Blink.Window()
body!(w, interface)
```
"""
function manualfit!(
    model::AbstractModel, data::NamedTuple{<:Any,<:Tuple{Vararg{<:AbstractVector}}};
    observations=[],
    throttle=0.1,
    kwargs...
)
    InteractModel(model; throttle=throttle) do updated_model
        @show params(updated_model)
        ModelParameters.setparent!(model, updated_model)
        manualfit(stripparams(updated_model), (observations, data); kwargs...)
    end
end

function manualfit(model, (observations, data); kwargs...)
    predictions = combinelayers(model, data)
    p = plot(first(data), predictions; label="predicted", legend=false, kwargs...)
    scatter!(p, observations; label="observed")
    return p
end

"""
    mapfit!(model::Model, modelkwargs;
            occurrence=[],
            precomputed=nothing,
            throttle=0.1,
            window=(Band(1),),
            kwargs...
    )

Fit a model to the map.

# Arguments

- `occurence`: a `Vector` of occurence locations, as `(lon, lat)` tuples.
- `modelkwargs`: are passed to the `mapgrowth` with the model.
- `mapkwargs`: are passed to the `plot` function the plots the `GeoArray`
- `throttle`: the response time of Interact.jl sliders.
- `window`: selects a window of the output to plot. By default this is `(Band(1),)`,
  which just removes the `Band` dimension from a `GeoArray`, if it exists.
- `kwargs`: passed to `Plots.scatter!`

# Example

```julia
model = Model(wiltstress, coldstress, heatstress)
interface = mapfit!(model, modelkwargs;
    occurrence=occurrence,
    precomputed=precomputed,
    throttle=0.2,
    window=(Band(1),),
    markershape=:cross,
    markercolor=:lightblue,
    markeropacity=0.4
)
display(interface)
```

To use the interface in a desktop app, use Blink.jl:

```julia; eval=false
using Blink
w = Blink.Window()
body!(w, interface)
```
"""
function mapfit!(model::AbstractModel, modelkwargs;
    mapkwargs=(),
    occurrence=[],
    precomputed=nothing,
    throttle=0.1,
    scatterkwargs...
)
    title = "GrowthMaps: mapfit interface"
    InteractModel(model; throttle=throttle, title=title) do updated_model
        ModelParameters.setparent!(model, updated_model)
        mapfit(
            stripparams(updated_model), (modelkwargs, occurrence, precomputed); 
            scatterkwargs...
        )
    end
end

function mapfit(model, (modelkwargs, occurrence, precomputed);
    window=(Band(1),),
    levels=10,
    mapkwargs=(),
    markercolor=:white,
    markersize=2.0,
    clims=(0.0, 0.25),
    scatterkwargs...
)
    output = mapgrowth(model; modelkwargs...)
    output = isnothing(precomputed) ? output : output .+ precomputed
    windowed = output[window...]
    p = plot(windowed; legend=:none, levels=levels, clims=clims, mapkwargs...)
    scatter(; t...) = scatter!(p, occurrence; t..., markercolor=markercolor, markersize=markersize, scatterkwargs...)
    if hasdim(windowed, Ti())
        for t in 1:length(dims(windowed, Ti()))
            scatter(; subplot=t)
        end
    else
        scatter()
    end
    return p
end
