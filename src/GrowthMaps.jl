module GrowthMaps

using Unitful, HDF5, Dates, GeoData

using GeoData: Time, rebuild
using Unitful: °C, K
using Base: tail

export growthrates

export AbstractRateModel, LowerStress, UpperStress, IntrinsicPopGrowth

include("models.jl")


"""
    growthrates(models, series::AbstractGeoSeries; 
                startdate=first(val(dims(series, Time))), nsubperiods=1, subperiod=Day(1), 
                period=Month(1), nperiods=12, constructor=identity)

Combine growth rates accross rate models and subperiods for all required periods.

The output is an array with the same dimensions as the passed in stack layers along a Time
dimension with a length of `nperiods`.
"""
function growthrates(models, series::AbstractGeoSeries;
                     startdate=first(val(dims(series, Time))),
                     nsubperiods=1,
                     subperiod=Day(1),
                     period=Month(1),
                     nperiods=12,
                     constructor=identity)
    # Allocate memory
    stack = series[Time<|At(startdate)]
    reqkeys = keys(models)
    # Make a memory backed stack using only the keys required for the model,
    # on the GPU if `constructor` is CuArray or similar
    layers = NamedTuple{reqkeys}(constructor.((stack[key] for key in reqkeys)))
    stackbuffer = rebuild(stack; parent=layers)
    mask = GeoData.mask(first(values(stackbuffer)))
    init = zero(first(values(stackbuffer)))

    # Setup output vector
    dates = val(dims(series, Time))
    periodstarts = [startdate + p * period for p in 1:nperiods]
    # Make a 3 dimensional GeoArray for output, adding the time dimension
    # to init (there should be a function for this in DimensionalData.jl - growdim?
    outputs = rebuild(init; parent=zeros(size(init)..., nperiods), dims=(dims(init)..., Time(periodstarts)))

    for p in 1:nperiods
        n = 0
        for sub in 1:nsubperiods
            # Subperiods are spaced evenly accross the period, with an offset so that a 
            # single subperiod falls in the middle of the period, and two are spaced 1/3 
            # and 2/3 along the period. They are rounded to start at eg. the start of a day
            # if typeof(subperiod) <: Day
            substart = round(periodstarts[p] + sub * (period ÷ (nsubperiods + 1)), subperiod)
            # Substract 1 second so Between won't include the first stack for the next subperiod
            subend = substart + subperiod - Second(1)
            subseries = series[Time<|Between(substart, subend)]
            length(subseries) > 0 || continue
            for t in size(subseries, Time)
                copy!(stackbuffer, subseries[t])
                outputs[p] .+= rate(models, stackbuffer)
                n += 1
            end
        end
        outputs[Time(p)] .*= mask ./ n
    end
    outputs
end

"""
Combine rates from all the model components for the current timestep data.

This recursive `@inline` method should fuse broadcasting accross all methods allowing 
the compiler to join them into a single CPU or GPU broadcast.
"""
@inline rate(models::Tuple, stack) = begin
    model = first(models)
    conditionalrate.(Ref(model), stack[keys(model)]) .+ rate(tail(models), stack)
end
@inline rate(models::Tuple{}, stack) = 0

@inline conditionalrate(model, val) = condition(model, val) ? rate(model, val) : 0

end # module
