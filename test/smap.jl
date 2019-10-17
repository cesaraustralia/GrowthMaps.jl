zipurl = "https://media.githubusercontent.com/media/cesaraustralia/GrowthMaps.jl/data/SMAP_aggregated27km.zip"
zipfilename = "SMAP_aggregated27km.zip"
folder = "SMAP_aggregated27km"

# Download SMAP zip
isfile(zipfilename) || download(zipurl)
# Unzip
run(`unzip -o $filename`)
filenames = readdir(folder)

# Separate wilting and surface temp filenames with regex
wilting = joinpath.(Ref(folder), filter(fn -> occursin(r"land_fraction_wilting", fn), filenames))
temp = joinpath.(Ref(folder), filter(fn -> occursin(r"surface_temp", fn), filenames))

# Get series dates with regex
df = DateFormat("yyyymmddTHHMMSS");
dates = DateTime.(replace.(temp, Ref(r".*_(\d+T\d+).tif" => s"\1")), Ref(df))

A = GDALarray(temp[1])
stacks = [GDALstack((land_fraction_wilting=wilting[i], surface_temp=temp[i]); 
                    dims=dims(A)) for i in 1:length(temp)]
stacks[1][:surface_temp]
series = GeoSeries(stacks, (Time(dates),));
series[Time(Near(DateTime(2016,1)))][:land_fraction_wilting] |> plot
stack = series[Near(DateTime(2016,1))][:surface_temp] |> plot

p = 3.377850e-01
ΔH_A = 3.574560e+04cal/mol
ΔH_L = -1.108990e+05cal/mol
ΔH_H = 3.276604e+05cal/mol
Thalf_L = 2.359187e+02K
Thalf_H = 2.991132e+02K
T_ref = K(25.0°C)
R = Unitful.R
growth = SchoolfieldIntrinsicGrowth(:surface_temp, p, ΔH_A, ΔH_L, 
                     ΔH_H, Thalf_L, Thalf_H, T_ref, R);
coldthresh = 7.0°C |> K  # Enriquez2017
coldmort = -log(1.00) * K^-1
coldstress = LowerStress(:surface_temp, coldthresh, coldmort)

heatthresh = 30.0°C |> K # Kimura2004
heatmort = -log(1.15) * K^-1
heatstress = UpperStress(:surface_temp, heatthresh, heatmort)

wiltthresh = 0.5 # default?
wiltmort = -log(1.1);
wiltstress = UpperStress(:land_fraction_wilting, wiltthresh, wiltmort);

model = growth, coldstress, heatstress, wiltstress

# TODO test output
output = mapgrowth(model, series;
                   startdate=DateTime(2016,1),
                   nperiods=12,
                   period=Month(1),
                   subperiod=Day(1))

using Plots, Statistics
max.(output[Time<|Near(DateTime(2016, 5))], 0) |> plot
