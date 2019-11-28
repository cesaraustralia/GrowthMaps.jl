"""
Passed to an optimiser to facilitate fitting any `RateModel`,
without boilerplate or methods rewrites.
"""
mutable struct ModelWrapper{M<:RateModel}
    model::M
end

"""
Update the model with passed in params, and run [`rate`] over the independent variables.
"""
(f::ModelWrapper)(xs, params) =
    GrowthMaps.conditionalrate.(Ref(reconstruct(f.model, params, Real)), xs)

"""
    fitlsq(model::RateModel, xs::AbstractArray, ys::AbstractArray)

Fit a model to data with least squares regression, using `curve_fit` from
LsqFit.jl. The passed in model should be initialised with sensible defaults,
these will be used as the initial parameters for the optimization.

Any (nested) `Real` fields on the struct are flattened to a parameter vector using
[Flatten.jl](http://github.com/rafaqz/Flatten.jl). Fields can be marked to ignore
using the `@flattenable` macro from [FieldMetadata.jl](http://github.com/rafaqz/FieldMetadata.jl).

## Arguments
- `model`: Any constructed [`RateModel`](@ref), including custom models.
- `xs`: AbstactArray of independent variables (model input values)
- `ys`: AbstactArray of dependent variables (model output values)

## Returns
The model with fitted parameters
"""
fit(model::ModelWrapper, xs::AbstractArray, ys::AbstractArray) = begin
    wrapper.model = fit(wrapper.model, xs, ys)
    wrapper
end
fit(model::RateModel, xs::AbstractArray, ys::AbstractArray) = begin
    fit = curve_fit(ModelWrapper(model), xs, ys, [flatten(model, Real)...])
    reconstruct(model, fit.param)
end


"""
Returns the wrapper holding the fitted model
"""
manualfit!(wrapper::ModelWrapper, xs, ys; throttlelen=0.1, plotrange=270.0K:0.1K:320K) = begin
    params = Flatten.flatten(val(wrapper), Number)
    fnames = fieldnameflatten(val(wrapper), Number)
    bounds = metaflatten(val(wrapper), FieldMetadata.bounds, Number)
    ranges = buildrange.(bounds, params)
    parents = parentnameflatten(val(wrapper), Number)
    descriptions = metaflatten(val(wrapper), FieldMetadata.description, Number)
    attributes = ((p, n, d) -> Dict(:title => "$p.$n: $d")).(parents, fnames, descriptions)

    plotobs = Observable(plotmodel(wrapper.model, xs, ys, params, plotrange))

    sliders = make_slider.(params, fnames, ranges, attributes)
    slider_obs = map((s...) -> s, throttle.(throttlelen, observe.(sliders))...)
    on(slider_obs) do params
        wrapper.model = Flatten.reconstruct(wrapper.model, params, Number)
        plotobs[] = plotmodel(wrapper.model, xs, ys, params, plotrange)
    end

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

    vbox(vbox(slider_groups...), plotobs)
end

plotmodel(model, xs, ys, params, plotrange) = begin
    predictions = map(x -> GrowthMaps.rate(model, x), plotrange)
    p = plot(plotrange, predictions; label="fitted", legend=false)
    scatter!(p, xs, ys; label="observed")
    p
end

make_slider(val, lab, rng, attr) = slider(rng; label=string(lab), value=val, attributes=attr)

buildrange(bounds::Tuple, val::T) where T =
    T(bounds[1]):(T(bounds[2])-T(bounds[1]))/1000:T(bounds[2])

electronfit(app; zoom=1.0) = begin
    ui = app(nothing)
    w = Window(Dict("webPreferences"=>Dict("zoomFactor"=>zoom)));
    # Blink.AtomShell.@dot w webContents.setZoomFactor($zoom)
    body!(w, ui);
end
