"""
    mapgrowth(model; series::AbstractGeoSeries, tspan::AbstractRange)
    mapgrowth(layers...; kwargs...)

Combine growth rates and stressors for all layers,
in all files in each `series` falling in each period of `tspan`.

## Arguments
- `model`: A `Model` or `Tuple` of `Layer` components. 
  `Layer`s can also be passed in as separate arguments.

## Keyword Arguments

- `series`: any AbstractGeoSeries from [GeoData.jl](http://github.com/rafaqz/GeoData.jl)
- `tspan`: `AbstractRange` for the timespan to run the layers for.
  This will become the `index` values of the output `GeoArray` time-dimension `Ti`.
- `arraytype`: An array constructor to apply to data once it's loaded from disk.
  The main use case for this is a `GPUArray` such as `CuArray` which will result in all
  computations happening on the GPU, if you have one.

Using multiple models in a `NamedTuple` can be an order of magnitude faster than
running models separately -- especially when `arraytype=CuArray` or similar.
In this configuration, all data required by all layers in all models will be loaded
from disk only once and copied to the GPU only once. This is useful for doing all kinds 
of sensitivity analysis and model comparison. 

##  Returns 

A `GeoArray` or `NamedTuple` of `GeoArray`, with the same dimensions as the passed-in 
stack layers, and an additional `Ti` (time) dimension matching `tspan`.
"""
mapgrowth(model::Model; kwargs...) = mapgrowth(parent(model); kwargs...)
mapgrowth(layers...; kwargs...) = mapgrowth(layers; kwargs...)
function mapgrowth(layers::Tuple; 
    series::AbstractGeoSeries, tspan::AbstractRange, arraytype=Array
)
    period = step(tspan); nperiods = length(tspan)
    startdate, enddate = first(tspan), last(tspan)
    required_keys = Tuple(union(map(l -> tuple(union(keys(l))), models)...)[1])

    # Copy only the required keys to a memory-backed stack
    stack = GeoStack(deepcopy(first(series)); keys=required_keys)

    A = first(values(stack));
    missingval = eltype(A)(NaN)

    # Replace false with NaN
    mask = map(x -> x ? eltype(A)(x) : missingval, boolmask(A)) |> parent |> arraytype
    stackbuffer = GeoData.modify(arraytype, stack)

    # Make a 3 dimensional GeoArray for output, adding the time dimension
    ti = Ti(tspan; mode=Sampled(Ordered(), Regular(period), Intervals(Start())))
    outdims = (dims(A)..., ti)
    outA = arraytype(zeros(eltype(A), size(A)..., nperiods))
    output = GeoArray(outA, outdims; name=:growthrate, missingval=missingval)

    runperiods!(outputs, stackbuffer, series, mask, models, tspan)

    # Return a GeoArray wrapping a regular Array, not arraytype
    map(o -> GeoData.modify(Array, o), outputs)
end

function runperiods!(outputs::NamedTuple, stackbuffer, series, mask, models::NamedTuple, tspan)
    period = step(tspan); nperiods = length(tspan)
    println("Running for $(1:nperiods)")
    for p in 1:nperiods
        n = 0
        periodstart = tspan[p]
        periodend = periodstart + period
        println("\n", "Processing period between: $periodstart and $periodend")

        # We don't use `Between` as it might unintentionally cut off the
        # last time if it partially extends beyond the period.
        # So we just work with time as Points using `Where`.
        subseries = series[Ti(Where(t -> t >= periodstart && t < periodend))]
        for t in 1:size(subseries, Ti)
            println("\n    ", val(dims(subseries, Ti))[t])
            # Copy the arrays we need from disk to the buffer stack
            copy!(stackbuffer, subseries[t])
            map(outputs, models, keys(outputs)) do output, model, key
                length(models) > 1 && println("        For $key")
                # For some reason now this is broken with DD getindex, view is a workaround
                parent(view(output, Ti(p))) .+= combinelayers(model, stackbuffer)
            end
            n += 1
        end
        if n > 0
            map(outputs) do output
                parent(view(output, Ti(p))) .*= mask ./ n
            end
        else
            @warn ("No files found for the $period period starting $periodstart")
            map(outputs) do output
                parent(view(output, Ti(p))) .*= mask
            end
        end
    end
end

@inline combinelayers(layer, stackbuffer) = combinelayers((layer,), stackbuffer)
@inline combinelayers(layers::Tuple, stackbuffer::AbstractGeoStack) =
    conditionalrate.(Ref(first(layers)), parent(stackbuffer[keys(first(layers))])) .+ combinelayers(tail(layers), stackbuffer)
@inline combinelayers(layers::Tuple{T}, stackbuffer::AbstractGeoStack) where T =
    conditionalrate.(Ref(first(layers)), parent(stackbuffer[keys(first(layers))]))
@inline combinelayers(layers::Tuple, vals::NamedTuple) =
    conditionalrate.(Ref(first(layers)), vals[keys(first(layers))]) .+ combinelayers(tail(layers), vals)
@inline combinelayers(layers::Tuple{T}, vals::NamedTuple) where T =
    conditionalrate.(Ref(first(layers)), vals[keys(first(layers))])


# Special-case Celcius: some data sets are in °C, but layers
# Can never use °C, so we convert it.
maybetoK(val) = unit(val) == u"°C" ? val |> K : val

@inline conditionalrate(l::Layer, val) = conditionalrate(model(l), maybetoK(unit(l) * val))
@inline conditionalrate(l::Layer, val::Unitful.AbstractQuantity) =
    conditionalrate(model(l), maybetoK(val))
@inline conditionalrate(model::RateModel, val) =
    condition(model, val) ? rate(model, val) : zero(rate(model, val))
