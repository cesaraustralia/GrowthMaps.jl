# GrowthMaps

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://cesaraustralia.github.io/GrowthMaps.jl/stable/)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://cesaraustralia.github.io/GrowthMaps.jl/dev/)
[![Build Status](https://travis-ci.org/cesaraustralia/GrowthMaps.jl.svg?branch=master)](https://travis-ci.org/cesaraustralia/GrowthMaps.jl)
[![codecov.io](http://codecov.io/github/cesaraustralia/GrowthMaps.jl/coverage.svg?branch=master)](http://codecov.io/github/cesaraustralia/GrowthMaps.jl?branch=master)

GrowthMaps.jl produces gridded growth rates from environmental data, by processing 
them with growth and stress models, following the method outlined in Maino et
al, _"Forecasting the potential distribution of the invasive vegetable leafminer
using ‘top-down’ and ‘bottom-up’ models"_ (in press). 

They are intended to be (and already practically used as) a replacement for CLIMEX and 
similar tools. Different from CLIMEX is that results arrays have units of growth/time. 
Another useful property of these models is that growth rate layers can be added and 
combined arbitrarily.

A primary use-case for GrowthMaps layers is in for calculating growth-rates for 
![Dispersal.jl](https://github.com/cesaraustralia/Dispersal.jl).

For data input, this package leverages
[`GeoData.jl`](http://github.com/rafaqz/GeoData.jl) to import stacks of
environmental data from many different sources, loaded lazily in sequence to
minimise memory use. These can also be loaded and run on a GPU.

## Example

```julia
using GrowthMaps, GeoData, HDF5, CUDA, Unitful

# Define a growth model
p = 3e-01
ΔH_A = 3e4cal/mol
ΔH_L = -1e5cal/mol
ΔH_H = 3e5cal/mol
Thalf_L = 2e2K
Thalf_H = 3e2K
T_ref = K(25.0°C)
growthmodel = SchoolfieldIntrinsicGrowth(p, ΔH_A, ΔH_L, Thalf_L, ΔH_H, Thalf_H, T_ref)

# Wrap the model with a data layer with the key
# to retreive the data from, and the Unitful.jl units.
growth = Layer(:surface_temp, K, growthmodel)
```

Now we will use GeoData.jl to load a series of [SMAP](https://smap.jpl.nasa.gov/) 
files lazily, and GrowthMaps.jl will load them to an Nvida GPU just in time for processing:

```julia

path = "your_SMAP_folder"
# Load 100s of HDF5 files lazyily with GeoData.jl
series = SMAPseries(path)
# Set the timespan you want layers for
tspan = DateTime(2016, 1):Month(1):DateTime(2016, 12)
# Use and Nvidia GPU for computations
arraytype = CuArray

# Run the model
output = mapgrowth(growth;
    series=aggseries,
    tspan=tspan,
    arraytype=arraytyps,
)

# Plot the first timestep
output[Ti(1)] |> plot
```

GrowthMaps.jl can run this growth model over thousands of HDF5 files in minutes, 
on a regular desktop with a GPU, although a CPU alone is not too much slower. 
You can also use an memory-backed arrays or NetCDF or GDAL files in `mapgrowth`.

The models can be chained together and run over multiple data layers simultaneously. 

See the [`Examples`](https://cesaraustralia.github.io/GrowthMaps.jl/dev/example/)
section in the documentation to get started. You can also work through the 
[example.jmd](https://github.com/cesaraustralia/GrowthMaps.jl/blob/master/docs/src/example.jmd) in atom
(with the language-weave plugin) or the
[notebook](https://github.com/cesaraustralia/GrowthMaps.jl/blob/gh-pages/dev/notebook/example.ipynb).


## Live Interfaces

GrowthMaps provides interfaces for manually fitting models where automated fits are not appropriate.

Model curves can be fitted to an `AbstractRange` of input data using `manualfit`:

![manualfit interface](https://github.com/cesaraustralia/GrowthMaps.jl/blob/media/manualfit.png?raw=true)

Observations can be fitted to a map. Aggregated maps `GeoSeries` can be used to fit models in real-time. 
See the examples for details.

![mapfit interface](https://github.com/cesaraustralia/GrowthMaps.jl/blob/media/mapfit.png?raw=true)
