
# struct CustomBulma <: InteractBase.WidgetTheme; end

# const custom_css = "/home/raf/julia/GrowthMaps/assets/interact.css"

# InteractBase.libraries(::CustomBulma) = [InteractBase.libraries(Interact.Bulma())...]
# InteractBase.registertheme!(:custombulma, CustomBulma())
# settheme!(CustomBulma())

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
- `obs`: `Vector` of observations, such as length 2 `Tuple`s of `Real`.
  Leave units off.

## Returns
An updated model with fitted parameters.
"""
fit(model, obs::AbstractArray) = begin
    fit = curve_fit(ModelWrapper(model), first.(obs), last.(obs), [flatten(model, Real)...])
    reconstruct(model, fit.param)
end

"""
    fit!(model::ModelWrapper, obs::AbstractArray)

Run `fit` on the contained model, and write the
updated model to the mutable modelwrapper.
"""
fit!(wrapper::ModelWrapper, obs::AbstractArray) = begin
    wrapper.model = fit(wrapper.model, obs)
    wrapper
end

"""
    manualfit!(wrapper::ModelWrapper, data::Array; obs=[],  kwargs...) =

Returns the wrapper holding the fitted model
"""
manualfit!(wrapper::ModelWrapper, ranges::NamedTuple{<:Any,<:Tuple{Vararg{<:AbstractVector}}};
           obs=[], kwargs...) =
    interface!(plotfit, wrapper, (obs, ranges); kwargs...)

plotfit(model, (obs, ranges)) = begin
    predictions = combinemodels(model, ranges)
    p = plot(first(ranges), predictions; label="predicted", legend=false)
    scatter!(p, obs; label="observed")
    p
end

"""
Fit a model to the map
"""
mapfit!(wrapper::ModelWrapper, modelkwargs; occurrence=[], precomputed=nothing, kwargs...) =
    interface!(plotmap, wrapper, (modelkwargs, occurrence, precomputed); kwargs...)

plotmap(model, (modelkwargs, occurrence, precomputed);
        window=(Band(1),), levels=10, markercolor=:white, markersize=2.0,
        clims=(0.0, 0.25), mapkwargs=(), scatterkwargs...) = begin
    output = mapgrowth(model; modelkwargs...)
    output = isnothing(precomputed) ? output : output .+ precomputed
    p = plot(output[window...]; legend=:none, levels=levels, clims=clims, mapkwargs...)
    for t in 1:length(dims(output, Ti()))
        scatter!(occurrence; subplot=t, markercolor=markercolor, markersize=markersize, scatterkwargs...)
    end
    p
end

interface!(f, wrapper::ModelWrapper, data;
           use=Number, ignore=Nothing, throttle=0.1, kwargs...) = begin
    plotobs = Observable(f(wrapper.model, data; kwargs...))
    sliders, slider_obs, slider_groups = build_sliders(wrapper.model, use, ignore, throttle)
    on(slider_obs) do params
        wrapper.model = Flatten.reconstruct(wrapper.model, params, use, ignore)
        plotobs[] = f(wrapper.model, data; kwargs...)
    end
    vbox(plotobs, vbox(slider_groups...))
end


build_sliders(model, use, ignore, _throttle) = begin
    params = Flatten.flatten(model, use, ignore)
    fnames = fieldnameflatten(model, use, ignore)
    bounds = metaflatten(model, FieldMetadata.bounds, use, ignore)
    ranges = buildrange.(bounds, params)
    parents = parentnameflatten(model, use, ignore)
    descriptions = metaflatten(model, FieldMetadata.description, use, ignore)
    attributes = ((p, n, d) -> Dict(:title => "$p.$n: $d")).(parents, fnames, descriptions)

    sliders = make_slider.(params, fnames, ranges, attributes)
    slider_obs = map((s...) -> s, throttle.(_throttle, observe.(sliders))...)

    group_title = nothing
    slider_groups = []
    group_items = []
    for i in 1:length(params)
        parent = parents[i]
        if group_title != parent
            group_title == nothing || push!(slider_groups, dom"div"(group_items...))
            group_items = Any[dom"h2"(string(parent))]
            group_title = parent
        end
        push!(group_items, sliders[i])
    end
    push!(slider_groups, dom"h2"(group_items...))

    sliders, slider_obs, slider_groups
end

make_slider(val, lab, rng, attr) =
    slider(rng; label=string(lab), value=val, attributes=attr)

buildrange(bounds::Tuple, val::T) where T =
    T(bounds[1]):(T(bounds[2])-T(bounds[1]))/1000:T(bounds[2])
