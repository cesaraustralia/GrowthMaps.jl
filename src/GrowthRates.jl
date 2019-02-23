module GrowthRates

using Unitful, HDF5, Dates
using Unitful: °C, K, d, h, ms
using Base: tail

export growthrates 

export AbstractRateModel, DryDays, HotDays, ColdDays, Wilting, IntrinsicPopGrowth

export AbstractDataSource, SMAP

include("models.jl")
include("smap.jl")


" Run all the stress models and collate to growth rate matrices"
function growthrates(models, data_path, start_date, end_date, timestep, sumspan; source=SMAP())
    # Get files and dates from filenames
    filepaths, dates = files_and_dates(source, data_path, start_date, end_date)

    # Allocate memory
    source_storage = allocate_source_storage(source, models, filepaths[1])
    init = source_storage[1]
    model_storage = allocate_model_storage(models, init)
    sub = arraysetup(similar(init, Bool))

    # Initialise counters
    day_counter = 0.0d
    spanstart = 0.0d
    spanend = spanstart + sumspan

    # Missings are slow. Convert them after calculations using a mask layer.
    mask = define_mask(source, init)
    
    # Setup output vector
    output = Array{Union{Missing,Float64},2}[]

    # Run for every date there is a file between the start and end dates 
    for i in 1:length(dates)
        println(dates[i])
        # We have allready initialised with data for the first date
        i > 1 && update_source_storage!(source, source_storage, models, filepaths[i])

        # Most days we use the whole day
        local frac = timestep


        # Get the number of days from the start date
        days = days_to_unit(dates[i] - dates[1])

        if days > spanend
            # Fractional day at the end of the span
            frac = days - spanend

            # Increment the number of days we have loaded by a fractional day
            day_counter += 1.0d - frac

            # Run models for last fractional day
            run_models!(model_storage, models, source_storage, 1d - frac)

            # Sum all stresses for the sumspan
            stresses = collate_stresses(model_storage, models, sumspan)

            # Remove units and mask output with missings
            push!(output, Array(stresses) .* mask)

            clear_model_storage!(model_storage)

            # Update the startpoint for the next sumspan
            spanstart = spanend 
            spanend = spanstart + sumspan
        else
            # Increment the number of days we have loaded by a full day
            day_counter += 1.0d
        end

        # Full day or fractional day at the beggining of the next summed timespan
        run_models!(model_storage, models, source_storage, frac)
    end
    output
end


"Allocate storage arrays for each submodel"
allocate_model_storage(models::Tuple, init) =
    (allocate_model_storage(models[1], init), allocate_model_storage(tail(models), init)...)
allocate_model_storage(::Tuple{}, init) = ()
allocate_model_storage(model, init) = 
    arraysetup([zero(typeof(1/scale(model))) for i in CartesianIndices(init)])

"Zero model arrays for the next timespan"
clear_model_storage!(model_storage::Tuple) = begin 
    model_storage[1] .= zero(eltype(model_storage[1])) 
    clear_model_storage!(tail(model_storage))
end
clear_model_storage!(model_storage::Tuple{}) = nothing

"Run all the models for the current timestep data, updating the storage array"
run_models!(model_storage::Tuple, models::Tuple, source_storage, frac) = begin
    model = models[1]
    source_data = source_storage[datakey(model)] * dataunit(model)
    model_storage[1] .+= subsettedrule.((model,), source_data) .* frac .* d^-1
    run_models!(tail(model_storage), tail(models), source_storage, frac)
end
run_models!(model_storage::Tuple{}, models::Tuple{}, source_storage, frac) = ()

@inline subsettedrule(model, val) = subset(model, val) ? rule(model, val) : val
    
"Combine all stresses for the timespan"
@inline collate_stresses(model_storage::Tuple, models::Tuple, sumspan) = 
    collate_stresses(model_storage[1], models[1], sumspan) .+ collate_stresses(tail(model_storage), tail(models), sumspan)
@inline collate_stresses(::Tuple{}, ::Tuple{}, sumspan) = 0
@inline collate_stresses(model_storage, model, sumspan) = 
    model_storage .* model.scale ./ sumspan .* d


"Convert julia Date days to Unitful days. We need floating point math."
days_to_unit(day::Day) = day.value * d
days_to_unit(millisecond::Millisecond) = millisecond.value * ms

""" 
Override this method to use a different array backend

To run on CUDA: 
```
arraysetup(a) = CuArrays(a)
```
"""
arraysetup(a) = a

end # module
