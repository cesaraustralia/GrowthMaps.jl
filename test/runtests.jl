using GrowthMaps, GeoData, HDF5, Dates, Unitful, Test
using GeoData: Time
using Unitful: °C, K, hr, d

month = Second(356.25*24*60*60/12) 
period = month
subperiod = Day(1)

# parameters
p25 = 3.377850e-01
H_A = 3.574560e+04K
H_L = -1.108990e+05K
T_0_5L = 2.359187e+02K
H_H = 3.276604e+05K
T_0_5H = 2.991132e+02K
tref = K(25.0°C)
R = 1.987u"kCal/K"

lowersm = 0.0 # no data
uppersm = 1.0 # no data
lowerct = 7.0°C |> K  # Enriquez2017
lowerctm = -log(1.00) * K^-1
upperct = 30.0°C |> K # Kimura2004
upperctm = -log(1.15) * K^-1
lowerwilt = 0.5 # default?
wiltingm = -log(1.1)

# models
growth = IntrinsicPopGrowth(:surface_temp, 1/K, p25, H_A, H_L, T_0_5L, H_H, T_0_5H, tref)
cold = LowerStress(:surface_temp, lowerct, lowerctm)
hot = UpperStress(:surface_temp, upperct, upperctm)
wilt = UpperStress(:land_fraction_wilting, lowerwilt, wiltingm)
models = (growth, cold, hot, wilt)

dims = Lon(1:4), Lat(1:5)
temp = GeoArray(ones(4, 5), dims)
wilt = GeoArray(ones(4, 5), dims)
stack = GeoStack(NamedTuple{(:surface_temp, :land_fraction_wilting)}((temp, wilt)), (),(),nothing)
series = GeoSeries([stack, stack], Time.((DateTime(2001), DateTime(2002))))
@which 
d = dims(series);
DimensionalData.dims(series::GeoSeries) = series.dims
series.dims

output = growthrates(models, series, period=period, subperiod=subperiod)



# Integration

struct TestModel{X,S,L} <: AbstractStressModel{X}
    scale::S
    lowersm::L
end
@inline rule(m::DryDays, x) = abs(x - m.lowersm)
@inline condition(m::DryDays, x) = x < m.lowersm





# Models
