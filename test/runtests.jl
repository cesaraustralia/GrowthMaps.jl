using GrowthMaps, GeoData, HDF5, Dates, Unitful, Test
using GeoData: Time, rebuild
using Unitful: °C, K, hr, d, mol, cal
using GrowthMaps: rate, condition, conditionalrate, period_startdates, subset_startdates
using Dates

# Models
dimz = Lat(10:20), Lon(100:130)

lowerdata = GeoArray([1 2 3
                      4 5 6], dimz)
key = :lower
threshold = 5K
mortalityrate = -1/K
lower = LowerStress(key, threshold, mortalityrate)
@test keys(lower) == :lower
@test condition.(Ref(lower), lowerdata) == [true true true
                                            true false  false]
@test conditionalrate.(Ref(lower), lowerdata) == [-4 -3 -2
                                                  -1  0  0]
upperdata = GeoArray([9 8 7
                      6 5 4], dimz)
key = :upper
threshold = 6K
mortalityrate = -2/K
upper = UpperStress(key, threshold, mortalityrate)
@test keys(upper) == :upper
@test condition.(Ref(upper), upperdata) == [true  true  true
                                            false false false]
@test conditionalrate.(Ref(upper), upperdata) == [-6 -4 -2
                                                   0  0  0]

key = :temp
p = 0.3
ΔH_A = 2e4cal * mol^-1
ΔH_L = -1e5cal * mol^-1
ΔH_H = 3e5cal * mol^-1
T_halfL = 250.0K
T_halfH = 300.0K
T_ref = K(25.0°C)
R = Unitful.R

tempdata = GeoArray([270.0 280.0 290.0
                     300.0 310.0 320.0], dimz)
growth = IntrinsicGrowth(key, p, ΔH_A, ΔH_L, ΔH_H, T_halfL, T_halfH, T_ref, R)
@test keys(growth) == :temp
rate.(Ref(growth), tempdata)
condition.(Ref(growth), tempdata)
conditionalrate.(Ref(growth), tempdata)

# Combined models

stack = GeoStack(NamedTuple{(:lower, :upper, :temp)}((lowerdata, upperdata, tempdata)))

@test conditionalrate((lower, upper), stack) == [-10 -7 -4
                                                  -1  0  0]
conditionalrate((growth, lower, upper), stack)

stack2 = GeoStack(NamedTuple{(:lower, :upper, :temp)}((lowerdata, upperdata, tempdata)))
stack3 = GeoStack(NamedTuple{(:lower, :upper, :temp)}((lowerdata, upperdata, tempdata)))
stack4 = GeoStack(NamedTuple{(:lower, :upper, :temp)}((lowerdata, upperdata, tempdata)))

dimz = (Time([DateTime(2016, 1, 1, 9), DateTime(2016, 1, 16, 15), 
              DateTime(2016, 2, 1, 10), DateTime(2016, 2, 16, 14)]),)
series = GeoSeries([stack, stack2, stack3, stack4], dimz)
model = (growth, lower, upper)

nperiods = 2
period = Month(1)
startdate = DateTime(2016)
enddate = startdate + period * nperiods

@test period_startdates(startdate, period, nperiods) == [DateTime(2016, 1, 1), DateTime(2016, 2, 1)] 

periodstart = DateTime(2016, 1)
periodend = DateTime(2016, 2)
substarts = vcat([[startdate + Month(m) + Day(d) for d in 0:10:30] for m in 0:1]...)
subs = subset_startdates(periodstart, periodend, substarts)
@test subs == DateTime.(2016, 1, [1, 11, 21, 31])
subs = subset_startdates(periodstart + Month(1), periodend + Month(1), substarts)
@test subs == DateTime.(2016, 2, [1, 11, 21])
substarts = periodstart:Day(1):periodend
subs = subset_startdates(periodstart, periodend, substarts)
@test length(subs) == 31
@test subs[20] == DateTime(2016, 1, 20)

subperiod_starts = DateTime.(2016, [1,1,2,2], [1, 16, 1, 16])
subs = subset_startdates(startdate, startdate+Month(1), subperiod_starts)
@test subs == DateTime.(2016, [1,1], [1, 16])
subs = subset_startdates(startdate+Month(1), startdate+Month(2), subperiod_starts)
@test subs == DateTime.(2016, [2,2], [1, 16])

output = mapgrowth(model, series; 
                   startdate=startdate,
                   nperiods=4,
                   period=Day(15),
                   subperiod=Day(1),
                   subperiod_starts=subperiod_starts,
                   constructor=identity)
@test typeof(dims(output)) <: Tuple{Lat,Lon,Time}
@test length(val(dims(output, Time))) == 4

output = mapgrowth(model, series; 
                   startdate=startdate,
                   nperiods=2,
                   period=Month(1),
                   subperiod=Day(2),
                   subperiod_starts=subperiod_starts,
                   constructor=identity)
@test length(val(dims(output, Time))) == 2


using GrowthMaps, GeoData, Dates, Plots, Unitful, Pkg
using GeoData: Time
# Load some Unitful.jl units
using Unitful: °C, K
basedir = Pkg.dir("GrowthMaps")

series = SMAPseries(joinpath(basedir, "test"); 
                    window=(Lon<|Between(-125, -75), Lat(150:480)));
series[Near(Date(2016))][:surface_temp] |> plot

# Define parameters
p = 3.377850e-01
ΔH_A = 3.574560e+04cal/mol
ΔH_L = -1.108990e+05cal/mol
ΔH_H = 3.276604e+05cal/mol
Thalf_L = 2.359187e+02K
Thalf_H = 2.991132e+02K
T_ref = K(25.0°C)
R = Unitful.R

coldthresh = 7.0°C |> K  # Enriquez2017
coldmort = -log(1.00) * K^-1
heatthresh = 30.0°C |> K # Kimura2004
heatmort = -log(1.15) * K^-1
wiltthresh = 0.5 # default?
wiltmort = -log(1.1);

# Define model components
growth = IntrinsicGrowth(:surface_temp, p, ΔH_A, ΔH_L, 
                         ΔH_H, Thalf_L, Thalf_H, T_ref, R)
coldstress = LowerStress(:surface_temp, coldthresh, coldmort)
heatstress = UpperStress(:surface_temp, heatthresh, heatmort)
wiltstress = UpperStress(:land_fraction_wilting, wiltthresh, wiltmort);

# Combine components into a model
model = (growth, coldstress, heatstress, wiltstress)
# Some alternative models:
model = (wiltstress,)
model = (growth, wiltstress)
model = (coldstress, heatstress)
model = (growth,)
model = (heatstress,)
model = (coldstress,)
model = (wiltstress,)

## Run your model
constructor=identity
output = mapgrowth(model, series; period=Month(1), nperiods=12, subperiod=Day(1))

pyplot()
output[Time(1)] |> plot
plot(series[2][:land_fraction_wilting])

using JLD2
path = joinpath(basedir, "docs/build")
mkpath(path)
rates = output
@save joinpath(path, "growthrates.jld2") rates 
