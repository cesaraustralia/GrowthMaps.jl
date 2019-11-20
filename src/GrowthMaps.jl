module GrowthMaps

# Use the README as the module docs
@doc read(joinpath(dirname(@__DIR__), "README.md"), String) GrowthMaps

using ConstructionBase, 
      Dates, 
      FieldMetadata,
      Flatten,
      GeoData, 
      HDF5, 
      LsqFit, 
      Unitful 

import FieldMetadata: @flattenable, flattenable
using GeoData: Time, rebuild
using Unitful: °C, K
using Base: tail

export mapgrowth, fit

export RateModel, GrowthModel, SchoolfieldIntrinsicGrowth, StressModel, LowerStress, UpperStress

include("models.jl")
include("fit.jl")


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
mapgrowth(model, series; kwargs...) = mapgrowth((model,), series; kwargs...)
mapgrowth(model::Tuple, series::AbstractGeoSeries;
          nperiods=1,
          period=Month(1),
          startdate=first(bounds(series, Time)),
          enddate=startdate + period * nperiods,
          subperiod=Day(1),
          subperiod_starts=startdate:subperiod:enddate,
          constructor=identity) = begin

    stack = first(series)
    reqkeys = Tuple(union(keys(model)))
    # Make a memory backed stack using only the keys required for the model,
    # on the GPU if `constructor` is CuArray or similar
    # TODO: make an operation that reqbuilds with a subset of keys in GeoData
    arraygen = (reconstructparent(stack[key], constructor) for key in reqkeys)
    stackbuffer = GeoStack(stack; data=NamedTuple{reqkeys}(Tuple(arraygen)))
    A = first(values(stackbuffer));
    mask = GeoData.boolmask(A)
    # Create an init array without refdims or a name
    init = GeoArray(A; data=zero(parent(A)), name=Symbol(""), refdims=());

    dates = val(dims(series, Time))
    periodstarts = period_startdates(startdate, period, nperiods)
    # Setup output vector
    # Make a 3 dimensional GeoArray for output, adding the time dimension
    # to init (there should be a function for this in DimensionalData.jl - growdim?
    output = rebuild(init; data=zeros(size(init)..., nperiods), 
                     dims=(dims(init)..., Time(periodstarts)))

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
                output[Time(p)] .+= conditionalrate(model, stackbuffer)
                n += 1
            end
        end
        if n > 0
            output[Time(p)] .*= parent(mask) ./ n
        else
            println("    No files found for this period")
            output[Time(p)] .*= parent(mask)
        end
        sleep(0)
    end

    output
end

reconstructparent(A, constructor) = GeoArray(A; data=constructor(parent(A)))

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

@inline conditionalrate(models::Tuple, stackbuffer) = 
    conditionalrate.(Ref(first(models)), parent(stackbuffer[keys(first(models))])) .+ conditionalrate(tail(models), stackbuffer)
@inline conditionalrate(models::Tuple{}, stackbuffer) = 0
@inline conditionalrate(model::RateModel, val) = condition(model, val) ? rate(model, val) : 0

end # module
