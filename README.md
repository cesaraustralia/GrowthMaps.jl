# GrowthMaps

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://cesaraustralia.github.io/GrowthMaps.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://cesaraustralia.github.io/GrowthMaps.jl/dev)
[![Build Status](https://travis-ci.org/cesaraustralia/GrowthMaps.jl.svg?branch=master)](https://travis-ci.org/cesaraustralia/GrowthMaps.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/dpf055yo50y21g1v?svg=true)](https://ci.appveyor.com/project/cesaraustralia/growthraates-jl)
[![codecov.io](http://codecov.io/github/cesaraustralia/GrowthMaps.jl/coverage.svg?branch=master)](http://codecov.io/github/cesaraustralia/GrowthMaps.jl?branch=master)
[![Coverage Status](https://coveralls.io/repos/cesaraustralia/GrowthMaps.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/cesaraustralia/GrowthMaps.jl?branch=master)

GrowthMaps.jl produces gridded growth rates from gridded environmental data and
process with growth and stress models, following the method outlined in Maino et
al, _"Forecasting the potential distribution of the invasive vegetable leafminer
using ‘top-down’ and ‘bottom-up’ models"_ (in press).

These models can be added and combined arbitrarily, and can run on a CPU or GPU.

For data input, this package leverages
[`GeoData.jl`](http://github.com/rafaqz/GeoData.jl) to import stacks of
environmental data from many different sources, loaded lazily in sequence to
minimise memory use. 

See the [`examples`](https://cesaraustralia.github.io/GrowthMaps.jl/dev/examples)
section in the documentation to get started.
