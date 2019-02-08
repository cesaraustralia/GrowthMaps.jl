module GrowthRates

using Unitful, HDF5, Dates
using Unitful: °C, K, hr, d
using Base: tail

export growthrates 

export AbstractRateModel, DryDays, HotDays, ColdDays, Wilting, IntrinsicPopGrowth

export AbstractDataSource, SMAP


# Fix unitful handling of missing
import Base.*
*(::Missing, ::T) where {T<:Unitful.Units} = missing
*(::T, ::Missing) where {T<:Unitful.Units} = missing


include("models.jl")
include("models.jl")


# Intermediate storage
struct RateModelStorage{A}
    accumulator::A
end

RateModelStorage(model, init) = begin
    typ = typeof(1/scale(model))
    a = arraysetup([zero(typ) for i in CartesianIndices(init)])
    RateModelStorage(a)
end


# Framework

" Run all the stress models and collate to growth rate matrices"
function growthrates(models, data_path, start_date, end_date, timespan; source=SMAP())

    filepaths, dates = daterange_filenames(source, data_path, start_date, end_date)
    filepaths = filepaths

    # Allocate memory
    source_storage = initialise_source_storage(source, models, filepaths[1])
    init = source_storage[1]
    storage = initialise_model_storage(models, init)
    sub = arraysetup(similar(init, Bool))

    day_counter = 0.0d
    timestep = 1.0d
    dayfrac = timestep / 1.0d
    

    for i in 1:length(dates)
        println(dates[i])
        i > 1 && update_source_storage!(source, source_storage, models, filepaths[i])
        run_models!(storage, models, source_storage, sub)
        day_counter += timestep
    end

    stresses = collate_stresses(models, storage) ./ day_counter ./ 1d
    for (i, s) in enumerate(storage)
        s.accumulator .= zero(1/scale(models[i])) 
    end

    mask = define_mask(source, init)
    Array(stresses) .* mask
end


initialise_model_storage(models::Tuple, init) =
    (RateModelStorage(models[1], init), initialise_model_storage(tail(models), init)...)
initialise_model_storage(::Tuple{}, init) = ()


run_models!(storage::Tuple, models::Tuple, source_storage, sub) = begin
    m = models[1]
    modeldata = source_storage[needsname(m)] * needsunit(m)
    sub .= subset.(Ref(m), modeldata)
    storage[1].accumulator .+= rule.((m,), modeldata)
    run_models!(tail(storage), tail(models), source_storage, sub)
end
run_models!(storage::Tuple{}, models::Tuple{}, source_storage, sub) = ()

    
@inline collate_stresses(models::Tuple, storage::Tuple) = 
    storage[1].accumulator .* models[1].scale .+ collate_stresses(tail(models), tail(storage))
@inline collate_stresses(::Tuple{}, ::Tuple{}) = 0


""" 
Override this method to use a different array backend

To run on CUDA: 
```
arraysetup(a) = CuArrays(a)
```
"""
arraysetup(a) = a

end # module
