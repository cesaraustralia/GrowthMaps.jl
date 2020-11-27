"""
    mapgrowth(model; series::AbstractGeoSeries, tspan::AbstractRange)
    mapgrowth(layers...; kwargs...)

Combine growth rates accross layers and subperiods for all required periods.

## Arguments

- `model`: A `Model`, NamedTuple of `Model`,  or `Tuple` of `Layer` components. 

`Layer`s can also be passed in as separate arguments. A `NamedTuple` of models
will all be run from the same data source, reducing load time.

## Keyword Arguments

- `series`: any AbstractGeoSeries from [GeoData.jl](http://github.com/rafaqz/GeoData.jl)
- `tspan`: `AbstractRange` for the timespan to run the layers for.
  This will become the `index` values of the output `GeoArray` time-dimension `Ti`.

The output is a GeoArray with the same dimensions as the passed in stack layers, and a Time
dimension with a length of `nperiods`.
"""
mapgrowth(layers::Layer...; kw...) = mapgrowth(layers; kw...)
mapgrowth(layers::Tuple{<:Layer,Vararg}; kw...) = mapgrowth(Model(layers); kw...)
mapgrowth(model::Model; kw...) = mapgrowth((model,); kw...)[1]
mapgrowth(m1::Model, ms::Model...; kw...) = mapgrowth((m1, ms...); kw...)
function mapgrowth(models::Union{Tuple{<:Model,Vararg},NamedTuple{<:Any,Tuple{<:Model,Vararg}}}; 
    series::AbstractGeoSeries, tspan::AbstractRange, arraytype=Array
)
    models = map(stripparams ∘ parent, models)
    period = step(tspan); nperiods = length(tspan)
    startdate, enddate = first(tspan), last(tspan)
    required_keys = Tuple(union(map(keys ∘ _astuple, models)...))

    # Copy only the required keys to a memory-backed stack
    stack = deepcopy(GeoStack(first(series); keys=required_keys))

    A = first(values(stack));
    missingval = eltype(A)(NaN)

    # Replace false with NaN
    mask = map(x -> x ? eltype(A)(x) : missingval, boolmask(A)) |> parent |> arraytype
    stackbuffer = GeoData.modify(arraytype, stack)

    # Make a 3 dimensional GeoArray for output, adding the time dimension
    ti = Ti(tspan; mode=Sampled(Ordered(), Regular(period), Intervals(Start())))
    outdims = (dims(A)..., ti)
    outA = arraytype(zeros(eltype(A), size(A)..., nperiods))
    outputs = map(models) do m
        GeoArray(deepcopy(outA), outdims; name=:growthrate, missingval=missingval)
    end

    runperiods!(outputs, stackbuffer, series, mask, models, tspan)

    # Return a GeoArray wrapping a regular Array, not arraytype
    map(o -> GeoData.modify(Array, o), outputs)
end

function runperiods!(outputs, stackbuffer, series, mask, models, tspan)
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
            println("    ", val(dims(subseries, Ti))[t])
            # Copy the arrays we need from disk to the buffer stack
            copy!(stackbuffer, subseries[t])
            map(outputs, models, keys(outputs)) do output, model, key
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

_astuple(x::Tuple) = x
_astuple(x) = (x,)
