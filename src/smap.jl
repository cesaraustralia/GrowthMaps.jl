
abstract type AbstractDataSource end


struct SMAP <: AbstractDataSource end

function define_mask(::SMAP, init)
    init = Array(init)
    mask = Union{Missing,Int}[init[i] == -9999.0 ? missing : 1 for i in CartesianIndices(init)]
end

initialise_source_storage(::SMAP, models, filepath) = begin
    geophysical = load_geophysical(filepath)
    required = tuple(union(needsname.(models))...)
    rasters = readkey.(Ref(geophysical), required)
    NamedTuple{required}(rasters)
end

update_source_storage!(::SMAP, source_storage, models, filepath) = begin
    geophysical = load_geophysical(filepath)
    for i = 1:length(source_storage) 
        source_storage[i] .= readkey(geophysical, keys(source_storage)[1])
    end
end

load_geophysical(filepath) = h5open(filepath)["Geophysical_Data"]
    
readkey(x, key) = arraysetup(read(x[string(key)]))

daterange_filenames(::SMAP, data_path, start_date, end_date) = begin
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
