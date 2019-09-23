module GrowthMaps
# Use the README as the module docs
@doc read(joinpath(dirname(@__DIR__), "README.md"), String) GrowthMaps

using Unitful, HDF5, Dates, GeoData

using GeoData: Time, rebuild
using Unitful: °C, K
using Base: tail

export mapgrowth

export AbstractRateModel, LowerStress, UpperStress, IntrinsicGrowth

include("models.jl")


"""
    growthrates(models, series::AbstractGeoSeries; 
                startdate=first(val(dims(series, Time))), nsubperiods=1, subperiod=Day(1), 
                period=Month(1), nperiods=12, constructor=identity)

Combine growth rates accross rate models and subperiods for all required periods.

The output is an array with the same dimensions as the passed in stack layers along a Time
dimension with a length of `nperiods`.
"""


function mapgrowth(model, series::AbstractGeoSeries;
                   nperiods=1,
                   period=Month(1),
                   startdate=first(val(dims(series, Time))),
                   enddate=startdate + period * nperiods,
                   subperiod=Day(1),
                   subperiod_starts=startdate:subperiod:enddate,
                   constructor=identity)
    # Allocate memory
    stack = first(series)
    reqkeys = keys(model)
    # Make a memory backed stack using only the keys required for the model,
    # on the GPU if `constructor` is CuArray or similar
    # TODO: make an operation that reqbuilds with a subset of keys in GeoData
    layers = NamedTuple{reqkeys}(constructor.((stack[key] for key in reqkeys)))
    stackbuffer = rebuild(stack; parent=layers)
    mask = GeoData.mask(first(values(stackbuffer)))
    init = zero(first(values(stackbuffer)))

    dates = val(dims(series, Time))
    periodstarts = period_startdates(startdate, period, nperiods)
    # Setup output vector
    # Make a 3 dimensional GeoArray for output, adding the time dimension
    # to init (there should be a function for this in DimensionalData.jl - growdim?
    outputs = rebuild(init; parent=zeros(size(init)..., nperiods), dims=(dims(init)..., Time(periodstarts)))

    for p in 1:nperiods
        periodstart = periodstarts[p]
        periodend = periodstart + period
        n = 0
        subs = subset_startdates(periodstart, periodend, subperiod_starts)
        for substart in subs
            subend = substart + subperiod - Second(1)
            subseries = series[Time<|Between(substart, subend)]
            # if size(subseries, Time) == 0 
                # println("Warning: no data found for the subperiod $substart to $subend")
            # end
            for t in 1:size(subseries, Time)
                copy!(stackbuffer, subseries[t])
                outputs[Time(p)] .+= conditionalrate(model, stackbuffer)
                n += 1
            end
        end
        outputs[Time(p)] .*= mask ./ n
    end
    outputs
end

period_startdates(startdate, period, nperiods) = [startdate + p * period for p in 0:nperiods-1]

subset_startdates(periodstart, periodend, substarts) = begin 
    # Get the first date in the period
    firstind = searchsortedfirst(substarts, periodstart)
    if firstind > length(substarts) || substarts[firstind] > periodend 
        error("No subperiods for period starting $periodstart")
    end
    lastind = searchsortedfirst(substarts, periodend) - 1
    substarts[firstind:lastind]
end

"""
Combine rates from all the model components for the current timestep data.

This recursive `@inline` method should fuse broadcasting accross all methods allowing 
the compiler to join them into a single CPU or GPU broadcast.
"""
@inline conditionalrate(models::Tuple, stack) = begin
    model = first(models)
    conditionalrate.(Ref(model), stack[keys(model)]) .+ conditionalrate(tail(models), stack)
end
@inline conditionalrate(models::Tuple{}, stack) = 0

@inline conditionalrate(model, val) = condition(model, val) ? rate(model, val) : 0

end # module
