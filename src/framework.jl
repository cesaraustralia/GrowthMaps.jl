"""
    mapgrowth(models, series::AbstractGeoSeries;
              startdate=first(val(dims(series, Time))), nsubperiods=1, subperiod=Day(1),
              period=Month(1), nperiods=12, constructor=identity)

Combine growth rates accross rate models and subperiods for all required periods.

## Arguments
- `models`: tuple of any RateModel compnents
- `series`: any AbstractGeoSeries from [GeoData.jl](http://github.com/rafaqz/GeoData.jl)

## Keyword Arguments
- `nperiods=12`: number of periods returned in the output
- `startdate=first(bounds(series, Time)))`: starting date of the sequence
- `enddate=startdate + period * nperiods`
- `period=Month(1)`: length of the period to output
- `subperiod=Day(1)`: length of the subperiods used to calculate output, such as a whole day to capture daily fluctuations.
- `constructor=identity`: Set to CuArray to process with CUDA on your GPU.

The output is a GeoArray with the same dimensions as the passed in stack layers, and a Time
dimension with a length of `nperiods`.
"""
mapgrowth(model, series; kwargs...) =
    mapgrowth((model,), series; kwargs...)
mapgrowth(wrapper::ModelWrapper, series; kwargs...) =
    mapgrowth(wrapper.model, series; kwargs...)
mapgrowth(model::Tuple, series::AbstractGeoSeries;
          nperiods=1,
          period=Month(1),
          startdate=first(GeoData.bounds(series, GeoData.Time)),
          enddate=startdate + period * nperiods,
          subperiod=Day(1),
          subperiod_starts=startdate:subperiod:enddate) = begin

    required_keys = Tuple(union(keys(model)))
    stackbuffer = GeoStack(first(series); keys=required_keys)
    A = first(values(stackbuffer));
    mask = map(x -> x ? eltype(A)(x) : eltype(A)(NaN), GeoData.boolmask(A))
    # Create an init array without refdims or a name
    init = GeoArray(A; data=zero(parent(A)), name="growthrate", refdims=())

    dates = val(dims(series, Time))
    periodstarts = period_startdates(startdate, period, nperiods)
    # Setup output vector
    # Make a 3 dimensional GeoArray for output, adding the time dimension
    # to init (there should be a function for this in DimensionalData.jl - growdim?
    output = GeoArray(init; data=zeros(size(init)..., nperiods),
                      dims=(dims(init)..., Time(periodstarts; grid=RegularGrid(; step=period))),
                      missingval=eltype(init)(NaN))

    println("Running for $(1:nperiods)")
    for p in 1:nperiods
        periodstart = periodstarts[p]
        println("\n", "Processing period starting: ", periodstart)
        periodend = periodstart + period
        n = 0
        subs = subset_startdates(periodstart, periodend, subperiod_starts)
        for substart in subs
            subend = min(substart + subperiod, periodend)
            subseries = series[Time<|Between(substart, subend - Second(1))]
            for t in 1:size(subseries, Time)
                println("    ", val(dims(subseries, Time))[t])
                copy!(stackbuffer, subseries[t])
                output[Time(p)] .+= combinemodels(model, stackbuffer)
                n += 1
            end
        end
        if n > 0
            output[Time(p)] .*= parent(mask) ./ n
        else
            println("    No files found for this period")
            output[Time(p)] .*= parent(mask)
        end
    end

    output
end

period_startdates(startdate, period, nperiods) = [startdate + p * period for p in 0:nperiods-1]

subset_startdates(periodstart, periodend, substarts) = begin
    # Get the first date in the period
    firstind = max(1, searchsortedfirst(substarts, periodstart))
    if firstind > length(substarts) || substarts[firstind] > periodend
        @warn "No subperiods found for period starting $periodstart"
    end
    local lastind = searchsortedlast(substarts, periodend)
    if substarts[lastind] < periodend
        substarts[firstind:lastind]
    else
        substarts[firstind:lastind-1]
    end
end

@inline combinemodels(models, stackbuffer) = combinemodels((models,), stackbuffer)
@inline combinemodels(models::Tuple, stackbuffer::AbstractGeoStack) =
    conditionalrate.(Ref(first(models)), parent(stackbuffer[keys(first(models))])) .+ combinemodels(tail(models), stackbuffer)
@inline combinemodels(models::Tuple{T}, stackbuffer::AbstractGeoStack) where T =
    conditionalrate.(Ref(first(models)), parent(stackbuffer[keys(first(models))]))
@inline combinemodels(models::Tuple, vals::NamedTuple) =
    conditionalrate.(Ref(first(models)), vals[keys(first(models))]) .+ combinemodels(tail(models), vals)
@inline combinemodels(models::Tuple{T}, vals::NamedTuple) where T =
    conditionalrate.(Ref(first(models)), vals[keys(first(models))])
@inline conditionalrate(model::Layer, val) = conditionalrate(model.model, val)
@inline conditionalrate(model::RateModel, val) = condition(model, val) ? rate(model, val) : zero(rate(model, val))
