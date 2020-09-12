# GrowthMaps

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://cesaraustralia.github.io/GrowthMaps.jl/stable/)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://cesaraustralia.github.io/GrowthMaps.jl/dev/)
[![Build Status](https://travis-ci.org/cesaraustralia/GrowthMaps.jl.svg?branch=master)](https://travis-ci.org/cesaraustralia/GrowthMaps.jl)
[![codecov.io](http://codecov.io/github/cesaraustralia/GrowthMaps.jl/coverage.svg?branch=master)](http://codecov.io/github/cesaraustralia/GrowthMaps.jl?branch=master)

GrowthMaps.jl produces gridded growth rates from environmental data, by processing
them with growth and stress models, following the method outlined in Maino et
al, _"Forecasting the potential distribution of the invasive vegetable leafminer
using ‘top-down’ and ‘bottom-up’ models"_
[(in press)](https://www.biorxiv.org/content/10.1101/866996v1).

They are intended to be (and already practically used as) a replacement for CLIMEX and
similar tools. Different from CLIMEX is that results arrays have units of growth/time.
Another useful property of these models is that growth rate layers can be added and
combined arbitrarily.

A primary use-case for GrowthMaps layers is in for calculating growth-rates for
![Dispersal.jl](https://github.com/cesaraustralia/Dispersal.jl).

For data input, this package leverages [`GeoData.jl`](http://github.com/rafaqz/GeoData.jl)
to import datasets from many different sources. Files are loaded lazily in sequence to
minimise memory, using the `GeoSeries` abstraction, that can hold "SMAP" HDF5 files,
NetCDFs, or `GeoStack`s of `tif` or other GDAL source files, or simply memory-backed
arrays. These data sources can be used interchangeably.

Computations can easily be run on GPUs, for example with `arraytype=CuArray`
to use CUDA for Nvidia GPUs.

## Example

Here we run a single growth model over SMAP data, on a Nvida GPU:

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
files lazily, and GrowthMaps.jl will load them to the GPU just in time for processing:

```julia
# Load 1000s of HDF5 files lazily using GeoData.jl
series = SMAPseries("your_SMAP_folder")

# Run the model
output = mapgrowth(growth;
    series=series,
    tspan=DateTime(2016, 1):Month(1):DateTime(2016, 12),
    arraytype=CuArray, # Use an Nvidia GPU for computations
)

# Plot every third month of 2016:
output[Ti(1:3:12)] |> plot
```

GrowthMaps.jl is fast.

As a rough benchmark, running a model using 3 3800*1600 layers from 3000 SMAP
files on a 7200rpm drive takes under 15 minutes on a desktop with a good GPU, 
like a GeForce 1080. The computation time is nearly irrelevant, running ten similar 
models takes essentially the same time as running one model. On a CPU, model run-time 
becomes more of a factor, but is still fast for a single model.

Additional improvements to the performance of `mapgrowth` may be found
largely by using solid-state and M.2. drives for data storage.

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
