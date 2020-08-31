# GrowthMaps

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://cesaraustralia.github.io/GrowthMaps.jl/stable/)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://cesaraustralia.github.io/GrowthMaps.jl/dev/)
[![Build Status](https://travis-ci.org/cesaraustralia/GrowthMaps.jl.svg?branch=master)](https://travis-ci.org/cesaraustralia/GrowthMaps.jl)
[![codecov.io](http://codecov.io/github/cesaraustralia/GrowthMaps.jl/coverage.svg?branch=master)](http://codecov.io/github/cesaraustralia/GrowthMaps.jl?branch=master)

GrowthMaps.jl produces gridded growth rates from gridded environmental data and
process with growth and stress models, following the method outlined in Maino et
al, _"Forecasting the potential distribution of the invasive vegetable leafminer
using ‘top-down’ and ‘bottom-up’ models"_ (in press). They are intended to 
be (and already practically used as) a replacement for CLIMEX and similar tools.
Different from CLIMEX is that results arrays have units of growth/time. 
Another useful property of these models is that growth rate layers can be added and 
combined arbitrarily.


![GrowthMaps output](https://github.com/cesaraustralia/GrowthMaps.jl/blob/gh-pages/dev/figures/example_19_1.png)

For data input, this package leverages
[`GeoData.jl`](http://github.com/rafaqz/GeoData.jl) to import stacks of
environmental data from many different sources, loaded lazily in sequence to
minimise memory use. These can also be loaded and run on a GPU.

See the [`Examples`](https://cesaraustralia.github.io/GrowthMaps.jl/dev/example/)
section in the documentation to get started.

## Live Interfaces

GrowthMaps provides interfaces for manually fitting models where automated fits are not appropriate.

Model curves can be fitted to an `AbstractRange` of input data using `manualfit`:

![manualfit interface](https://github.com/cesaraustralia/GrowthMaps.jl/blob/media/manualfit.png?raw=true)

Observations can be fitted to a map. Aggregated maps `GeoSeries` can be used to fit models in real-time. 
See the examples for details.

![mapfit interface](https://github.com/cesaraustralia/GrowthMaps.jl/blob/media/mapfit.png?raw=true)
