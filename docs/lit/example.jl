# # Modelling establishment potetial
#
# In this example we will calculate the expected population growth rates of
# Spotted Wing Drosophila (SWD) _D. suzukii_, for each month of the year on a 9km grid
# accross North America.
#
# We'll experiment with running the model with a number of different datasets.
# GrowthMaps.jl, via GeoData.jl faciliates using a wide range of input sources.
# First we will use aggregated data in tiff files, then aggregate them to in-memory
# files to run the inferface. Finally if you can download it, we'll use the
# (huge) SMAP dataset, which is in a custom HDF5 format.
#
# If you are using this as a .jmd file, it's best in
# to use atom with the "uber-juno" and "language-weave" plugins.
#
# ## Load some required packages
#
# These packages take care of loading and plotting data, and handling sci units and dates.
using GeoData, GrowthMaps, Plots, Unitful, UnitfulRecipes, Dates, CSV, DataFrames
basedir = realpath(joinpath(dirname(pathof(GrowthMaps)), "../docs"))

# ## Define model components

# First we'll define the growth model using `SchoolfieldIntrinsicGrowth`,
# based on Schoolfield (1981).

# ### Growth Model
#
# When defining model components, the first parameter is a `:symbol` for the
# required raster layer in the source data.
p = Param(3e-1; bounds=(3e-2, 3e0))
ΔH_A = Param(3e4; units=u"cal/mol", bounds=(3e3, 3e5))
ΔH_L = Param(-1e5; units=u"cal/mol", bounds=(-1e6, -1e4))
ΔH_H = Param(3e5; units=u"cal/mol", bounds=(3e4, 3e6))
Thalf_L = Param(2e2; units=u"K", bounds=(2e1, 2e3))
Thalf_H = Param(3e2; units=u"K", bounds=(3e1, 3e3))
T_ref = u"K"(25.0u"°C")

# Use the parameters in the SchoolfieldIntrinsicGrowth model
growthmodel = SchoolfieldIntrinsicGrowth(p, ΔH_A, ΔH_L, Thalf_L, ΔH_H, Thalf_H, T_ref)

# Now connect the model with a layer from the dataset, and give the layer units
growth = Layer(:surface_temp, u"K", growthmodel)

# Load lab observations
#
# If these are only estimated parameters, we can also fit the model to a a dataset
# of growth rate and temperature. 
#
# First extract our independent and dependent variables from the example CSV:
obsdata = CSV.File(joinpath(basedir, "swd_ecophys_data.csv"); 
    select=[:x_value, :y_value, :y_key]) |>
    DataFrame |>
    df -> filter(d -> d.y_key == "r_m", df) |>
    dropmissing

# Then extract the required data colummns, and convert temperature values from
# unitless Celcius to explit Kelvins, using Unitful.jl:
obsrate = obsdata.y_value
obstemp = obsdata.x_value .* u"°C" .|> u"K"
obs = collect(zip(obstemp, obsrate))

# Now we can fit the model. The `fitrate` function provides an easy way to fit the
# models in GrowthMaps or your own custom `RateModel`s, using the LsqFit.jl
# package:
model = Model(growth)
GrowthMaps.fit!(model, obs)

# Plot the fit against the data
# 
# Define a custom temperature range
temprange = (270.0:0.1:310.0)u"K"
# Plot the rate
temp_plot = plot(temprange; label="fitted") do x
    GrowthMaps.rate(stripparams(model), x) 
end
scatter!(temp_plot, obs; label="observed")

# Manual tweaking
#
# We can also try tweaking the fitting the model manually in a user interface.
# Model components are immutable (for performance reasons), so we wrap the model
# in a mutable wraper so we can use the results.
# We parametrise the model over the same temperature range that we are plotting,
# using the :surface_temp key that the model requires:
tempdata = (surface_temp=temprange,)
fit_interface = manualfit!(model, tempdata; obs=obs)

# To use the it in a desktop app, use Blink.jl:
## using Blink
## w = Blink.Window()
## body!(w, fit_interface)

# Note that `manualfit!` will also work with multiple model components
# that use the same source data, like `Model(growth, heatstress, coldstress)`.
#
# ## Load spatial data

# Later we can use real SMAP datasets using GeoData.jl SMAPseries loader.
# But downloading the dataset takes too long for an example.
# Instead we will download and unzip some lower resolution monthly
# data to use in the model:
dataurl = "https://media.githubusercontent.com/media/cesaraustralia/GrowthMaps.jl/data/SMAP_aggregated27km.zip"
zipfilepath = joinpath(basedir, "SMAP_aggregated27km.zip")
unzippedfolder = joinpath(basedir, "SMAP_aggregated27km")
isfile(zipfilepath) || download(dataurl, zipfilepath)
run(`unzip -o $zipfilepath -d $basedir`);


# Get the paths for to all the wilting and surface-temp files
filenames = readdir(unzippedfolder)
wilting_filenames = filter(fn -> occursin(r"land_fraction_wilting", fn), filenames)
surface_temp_filenames = filter(fn -> occursin(r"surface_temp", fn), filenames)
wilting_paths = joinpath.(Ref(unzippedfolder), wilting_filenames)
surface_temp_paths = map(p -> joinpath(unzippedfolder, p), surface_temp_filenames)

# Get the dates the `surface_temp` files using regex.
# Then use them to create the time dimension
date_format = DateFormat("yyyymmddTHHMMSS");
dates = DateTime.(replace.(surface_temp_paths, Ref(r".*_(\d+T\d+).tif" => s"\1")), Ref(date_format))
timedim = Ti(dates; mode=Sampled(span=Regular(Hour(3))))

# We know the `land_fraction_wilting` files are for the same dates.
#
# Now we have the files and date seies, we can put together a series of
# GeoData.jl stacks to load lazily from disk while running `mapgrowth`.
#
# Make a NamedTuple of Vectors of String paths to the files
layers = (land_fraction_wilting=wilting_paths, surface_temp=surface_temp_paths)

# Create the series
tiffseries = GeoSeries(layers, (timedim,); mappedcrs=EPSG(4326))
tiffseries = map(x -> view(x, Band(1)), tiffseries)

# We can plot a layer from a file at some date in the series:
tiffseries[Ti(Near(DateTime(2016, 2)))] |> plot

# Set the time period to a month, and set the length of the subsample period to
# the available times over one day:
#
# ## Run a model over the spatial data
#
# Define a timespan range to run the model over:
tspan = DateTime(2016, 1):Month(1):DateTime(2016, 12)

# Then to start, we'll run a simple model that only calculates growth rate.
growth_output = mapgrowth(model;
    series=tiffseries, tspan=tspan,
)

# Then plot the results:
growth_output[Ti(1)] |> plot

# It doesn't really capture realistic population growth: there are no growth rates
# below zero. We need to add some stress models. Stress models model processes
# that produce negative deathrates in the population, as oposed to rate models,
# where the minimum growth rate is zero.
#
# Stress models have a threshold in K and mortality
# rate per degree K. `LowerStress` models stress below the given threshold, while
# `UpperStress` models stress induced above a threshold.
#
# We will define stress models for cold, heat and wilt stress.
# As SWD live on and around plants, we use the proportion of plants
# wilting as an indicater of stress induced by lack of moisture.
coldthresh = Param(ustrip(u"K", -10.0f0u"°C"); units=u"K", bounds=(240, 290))  # Stephens2015
coldmort = Param(-log(1.23f0); units=u"K"^-1, bounds=(0, 0.4)) # Stephens2015
coldstress = Layer(:surface_temp, u"K", LowerStress(coldthresh, coldmort))
log(1.23f0)

heatthresh = Param(ustrip(u"K"(30.0u"°C")); units=u"K", bounds=(280, 330)) # Kimura2004
heatmort = Param(-log(1.15); units= u"K"^-1, bounds=(0, 0.4))
heatstress = Layer(:surface_temp, u"K", UpperStress(heatthresh, heatmort))

wiltthresh = Param(0.5; bounds=(0, 1)) # default?
wiltmort = Param(-log(1.1); bounds=(0, 0.4));
wiltstress = Layer(:land_fraction_wilting, UpperStress(wiltthresh, wiltmort));

# Missing the aggregated data for this!
# wetthresh = 0.8f0
# wetmort = -10.0f0
# wetstress = Layer(:sm_surface_wetness, WetStress(wetthresh, wetmort));

# ## Analyse the stressors
#
# ## Heat stress
heat_output = mapgrowth(heatstress;
    series=tiffseries, tspan=tspan,
)
plot(heat_output[Ti(1)])

# ## Cold stress
cold_output = mapgrowth(coldstress;
    series=tiffseries, tspan=tspan,
)
plot(cold_output[Ti(1)])

# ## Wilting stress
wilt_output = mapgrowth(wiltstress;
    series=tiffseries, tspan=tspan,
)
plot(wilt_output[Ti(1)])

# ## Combine the stressors
stress_output = mapgrowth(coldstress, heatstress, wiltstress;
    series=tiffseries, tspan=tspan,
)
plot(stress_output[Ti(1)])

plot(stress_output[Ti(1)]; clims=(-1, 0))

# To build a more complex model, we can chain components together
# in a tuple, and again, run and plot them:
model_output = mapgrowth(growth, coldstress, heatstress, wiltstress;
    series=tiffseries, tspan=tspan,
)

# ## Show just rates above zero
#
# For January:
plot(model_output[Ti(Jan)]; color=:oleron, clims=(-0.3, 0.3))

# For July:
plot(model_output[Ti(July)]; color=:oleron, clims=(-0.3, 0.3))

# ## Take the mean growth rate
using Statistics
mean_growth = mean(model_output; dims=Ti)
plot(mean_growth; color=:oleron, clims=(-0.3, 0.3))

# ## Compare with observation data
#
# To compare out simulation with observations data,
# we'll first load them from a CSV file:
csvurl = "https://raw.githubusercontent.com/cesaraustralia/GrowthMaps.jl/data/Oersted_occurrence.csv"
csvfilename = joinpath(basedir, "Oersted_occurrence.csv")
isfile(csvfilename) || download(csvurl, csvfilename)

obs = CSV.File(csvfilename)
occurrence = collect(zip(obs.Longitude, obs.Latitude))

# ### Plot with a split colorscheme
#
# This means positive and negative growth-rates are clear
growth_plot = plot(mean_growth; color=:oleron, clims=(-0.3, 0.3))  
# Then scatter them on a map:
scatter!(growth_plot, occurrence; 
     label="obs", legend=:none, fillalpha=0,
     markeralpha=0.4,
     markersize=2.0, markercolor=:red, markershape=:cross,
)
savefig("growth_and_obs.png")

# ## Parametrising models using interactive maps

# If you need to adjust the model based on the distribution, this
# can be done live in the interface, as with the manual fit.
#
# But parametrising maps on large datasets is processor intensive, which inhibits
# interactive fedback. To reduce processing, we can aggregate the spatial data to
# a more manageable size.
#
# You can experiment with the `agg` size to compromise between quality and render time.
# Large values will look pixelated but will run fast.
#
# `Center()` simply takes the central cell. `Statistics.mean` would take the mean value.
agg = 8
aggseries = GeoData.aggregate(Center(), read(tiffseries), (X(agg), Y(agg)))

# Then fit the other components. `throttle` will dictate how fast the interface
# updates. Make it larger on a slow machine, smaller on a faster one.
# The `window` argument lets us select a window of the output to view, using
# DimensionalData.jl/GeoData.jl `Dimension`s and `Selector`s. Here we only plot
# the first time: but you can select any/multiple, and select them with
# `Near(DateTime(2016,2))` or similar.
whole_year = DateTime(2016, 1):Year(1):DateTime(2016, 12)
stress_model = Model((coldstress, heatstress, wiltstress))
mapfit_interface = mapfit!(model, (series=aggseries, tspan=whole_year);
    occurrence=occurrence,
    throttle=0.1,
    colorbar=true,
    clims=(-2, 0),
    markershape=:cross,
    markercolor=:red,
    markeropacity=0.4,
    window=(Ti(1),),
)

# In atom, this will show the interface in the plot pane:
display(mapfit_interface)

# To use the it in a desktop app, use Blink.jl:
## using Blink
## w = Blink.Window()
## body!(w, mapfit_interface)

# And get the updated model components from the model:
# wiltstress, coldstress, heatstress = parent(stress_model)

# Now we will put together decent population growth maps
# using the higher resolutions data, and a monthly tiestep:
final_output = mapgrowth(growth, wiltstress, coldstress, heatstress;
    series=tiffseries,
    tspan=tspan,
);
plot(final_output[Ti(1:3:12)]; axis=false)


# ## Save the results
#
# Write the output as a NetCDF file
write("growthrates.nc", final_output)

# We can load it again with:
A = GeoArray("growthrates.nc")
plot(A[Ti(1)])


# Use the real SMAP dataset

# Using aggregated data loses some of the impact of the stress
# models, which respond to extremes, not averages. If you need to use
# the model for a serious application, run it on the real SMAP dataset.
#
# Unfortunately we can't just download the data for you with a script,
# as you need to log in to an account to have access.
#
# But GeoData.jl has a `SMAPseries()` constructor that will automate the
# whole process of loading SMAP HDF5 files from a folder once you have
# the data.
#
# We can also reduce the numer of days processed to 8 each month (or whatever)
# using the `Where` selector from DimensionalData.jl.

## using HDF5

# Set this to your own SMAP data folder
## smappath = "/home/raf/Storage/Data/SMAP/SMAP_L4_SM_gph_v4"

# Choose the dayes of month to use
## days = (1, 4, 8, 12, 16, 20, 23, 27) # 8 days per month

#days = 1:31 # all days
# Use the `Where` selector to choose dates in the time index by day of month
## smapseries = SMAPseries(smappath)[Ti(Where(d -> dayofmonth(d) in days))]

# Run the model:
## @time output = mapgrowth(model, wiltstress, coldstress, heatstress;
##    series=smapseries, tspan=tspan,
## )
## plot(output[Ti(1:3:12)]; axis=false)

# And save a netcdf:
## write(joinpath(basedir, "growthrates.nc"), output)
