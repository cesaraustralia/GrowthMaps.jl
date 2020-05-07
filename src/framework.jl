"""
    mapgrowth(models, series::AbstractGeoSeries;
              startdate=first(val(dims(series, Ti))), nsubperiods=1, subperiod=Day(1),
              period=Month(1), nperiods=12, constructor=identity)

Combine growth rates accross rate models and subperiods for all required periods.

## Arguments
- `models`: tuple of any RateModel compnents
- `series`: any AbstractGeoSeries from [GeoData.jl](http://github.com/rafaqz/GeoData.jl)

## Keyword Arguments
- `nperiods=12`: number of periods returned in the output
- `startdate=first(bounds(series, Ti)))`: starting date of the sequence
- `enddate=startdate + period * nperiods`
- `period=Month(1)`: length of the period to output
- `subperiod=Day(1)`: length of the subperiods used to calculate output, such as a whole day to capture daily fluctuations.
- `constructor=identity`: Set to CuArray to process with CUDA on your GPU.

The output is a GeoArray with the same dimensions as the passed in stack layers, and a Time
dimension with a length of `nperiods`.
"""
mapgrowth(; model, kwargs...) =
    mapgrowth(model; kwargs...)
mapgrowth(wrapper::ModelWrapper; kwargs...) =
    mapgrowth(wrapper.model; kwargs...)
mapgrowth(model...; kwargs...) =
    mapgrowth(model; kwargs...)
mapgrowth(model::Tuple; 
          series,
          period=Month(1),
          nperiods=1,
          startdate=first(bounds(series, Ti))) = begin
    println("main")
    enddate = startdate + period * nperiods
    required_keys = Tuple(union(keys(model)))
    # Copy only the required keys to a memory-backed stack
    stackbuffer = GeoStack(deepcopy(first(series)); keys=required_keys)
    A = first(values(stackbuffer));
    mask = map(x -> x ? eltype(A)(x) : eltype(A)(NaN), boolmask(A))

    periodstarts = periodstartdates(startdate, period, nperiods)
    # Setup output vector
    # Make a 3 dimensional GeoArray for output, adding the time dimension
    # to init (there should be a function for this in DimensionalData.jl - growdim?
    ti = Ti(periodstarts; mode=Sampled(Ordered(), Regular(period), Intervals(Start())))
    output = GeoArray(zeros(size(A)..., nperiods), (dims(A)..., ti);
        refdims=(),
        name="growthrate",
        missingval=eltype(A)(NaN)
    )

    println("Running for $(1:nperiods)")
    for p in 1:nperiods
        n = 0
        periodstart = periodstarts[p]
        periodend = periodstart + period
        println("\n", "Processing period between: $periodstart and $periodend")
        # We don't use `Between` as it might unintentionally cut off the
        # last time if it partially extends beyond the period.
        # So we jsut work with time as Points using `Where`.
        subseries = series[Ti(Where(t -> t >= periodstart && t < periodend))]
        for t in 1:size(subseries, Ti)
            println("    ", val(dims(subseries, Ti))[t])
            # Copy the arrays we need from disk to the buffer stack
            copy!(stackbuffer, subseries[t])
            output[Ti(p)] .+= combinemodels(model, stackbuffer)
            n += 1
        end
        if n > 0
            output[Ti(p)] .*= parent(mask) ./ n
        else
            @warn ("No files found for the $period period starting $periodstart")
            output[Ti(p)] .*= parent(mask)
        end
    end

    output
end

periodstartdates(startdate, period, nperiods) =
    [startdate + p * period for p in 0:nperiods-1]

@inline combinemodels(models, stackbuffer) = combinemodels((models,), stackbuffer)
@inline combinemodels(models::Tuple, stackbuffer::AbstractGeoStack) =
    conditionalrate.(Ref(first(models)), parent(stackbuffer[keys(first(models))])) .+ combinemodels(tail(models), stackbuffer)
@inline combinemodels(models::Tuple{T}, stackbuffer::AbstractGeoStack) where T =
    conditionalrate.(Ref(first(models)), parent(stackbuffer[keys(first(models))]))
@inline combinemodels(models::Tuple, vals::NamedTuple) =
    conditionalrate.(Ref(first(models)), vals[keys(first(models))]) .+ combinemodels(tail(models), vals)
@inline combinemodels(models::Tuple{T}, vals::NamedTuple) where T =
    conditionalrate.(Ref(first(models)), vals[keys(first(models))])
@inline conditionalrate(l::Layer, val) = 
    conditionalrate(model(l), unit(l) * val)
@inline conditionalrate(l::Layer, val::Quantity) = 
    conditionalrate(model(l), val)
@inline conditionalrate(model::RateModel, val) = 
    condition(model, val) ? rate(model, val) : zero(rate(model, val))

