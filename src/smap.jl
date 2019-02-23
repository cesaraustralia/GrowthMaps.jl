# SMAP data backenj

abstract type AbstractDataSource end


struct SMAP <: AbstractDataSource end

function define_mask(::SMAP, init)
    init = Array(init)
    mask = Union{Missing,Int}[init[i] == -9999.0 ? missing : 1 for i in CartesianIndices(init)]
end

allocate_source_storage(::SMAP, models, filepath) = begin
    # Get all the data keys reuired for all the models
    requiredkeys = tuple(union(datakey.(models))...)
    # Get all the data matrching the keys
    rasters = h5open(filepath) do data
        readkey.(Ref(data["Geophysical_Data"]), requiredkeys)
    end
    # Put data in a named tuple with the same keys
    NamedTuple{requiredkeys}(rasters)
end

update_source_storage!(::SMAP, source_storage, models, filepath) = begin
    requiredkeys = keys(source_storage)
    # Update arrays for each required key with the new dataset
    rasters = h5open(filepath) do data
        for i = 1:length(source_storage) 
            source_storage[i] .= readkey(data["Geophysical_Data"], requiredkeys[i])
        end
    end
end

files_and_dates(::SMAP, data_path, start_date, end_date) = begin
    filenames = readdir(data_path)
    pattern = Regex("_Vv401(1|0)_001.h5.iso.xml")
    datanames = [replace(f, "_Vv4011_001.h5.iso.xml" => "") for f in filenames if occursin(pattern, f)]
    smapformat = dateformat"\S\MAP_L4_\S\M_gph_yyyymmddTHHMMSS"
    dates = DateTime.(datanames, smapformat)
    # Find subset of all available dates between start and end dates
    subind = (dates .>= start_date) .& (dates .<= end_date)
    # Subset dates and filenames
    dates = dates[subind]
    filepaths = data_path .* datanames[subind] .* "_Vv4011_001.h5"
    filepaths, dates
end
    
readkey(x, key) = arraysetup(read(x[string(key)]))
