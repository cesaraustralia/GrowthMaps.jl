using GrowthMaps, GeoData, HDF5, Dates, Unitful, Test
using GeoData: Time, rebuild
using Unitful: °C, K, hr, d, mol, cal
using GrowthMaps: rate, condition, conditionalrate, period_startdates, subset_startdates
using Dates

dimz = Lat(10:20), Lon(100:130)

# Models
lowerdata = GeoArray([1 2 3
                      4 5 6], dimz)
upperdata = GeoArray([9 8 7
                      6 5 4], dimz)
tempdata = GeoArray([270.0 280.0 290.0
                     300.0 310.0 320.0], dimz)

key = :lower
threshold = 5K
mortalityrate = -1/K
lower = LowerStress(key, threshold, mortalityrate)
    
@testset "Lower stress" begin
    @test keys(lower) == :lower
    @test condition.(Ref(lower), lowerdata) == [true true true
                                                true false  false]
    @test conditionalrate.(Ref(lower), lowerdata) == [-4 -3 -2
                                                      -1  0  0]
end

key = :upper
threshold = 6K
mortalityrate = -2/K
upper = UpperStress(key, threshold, mortalityrate)

@testset "Upper stress" begin
    @test keys(upper) == :upper
    @test condition.(Ref(upper), upperdata) == [true  true  true
                                                false false false]
    @test conditionalrate.(Ref(upper), upperdata) == [-6 -4 -2
                                                       0  0  0]
end

key = :temp
p = 0.3
ΔH_A = 2e4cal * mol^-1
ΔH_L = -1e5cal * mol^-1
ΔH_H = 3e5cal * mol^-1
T_halfL = 250.0K
T_halfH = 300.0K
T_ref = K(25.0°C)
R = Unitful.R
growth = SchoolfieldIntrinsicGrowth(key, p, ΔH_A, ΔH_L, ΔH_H, T_halfL, T_halfH, T_ref, R)

@testset "Schoolfield Intrinsic growth" begin
    @test keys(growth) == :temp
    @test_broken condition.(Ref(growth), tempdata) # == TODO
    @test_broken conditionalrate.(Ref(growth), tempdata) # == TODO
end


@testset "Combined models" begin
    stack = GeoStack(NamedTuple{(:lower, :upper, :temp)}((lowerdata, upperdata, tempdata)))
    @test conditionalrate((lower, upper), stack) == [-10 -7 -4
                                                      -1  0  0]
    conditionalrate((growth, lower, upper), stack)
end

@testset "selection of dates for periods and subperiods" begin
    nperiods = 2
    period = Month(1)
    startdate = DateTime(2016)
    @testset "Test generation of start dates" begin
        @test period_startdates(startdate, period, nperiods) == [DateTime(2016, 1, 1), DateTime(2016, 2, 1)]
    end

    @testset "Test selection of start dates" begin
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
    end
end

@testset "Integration" begin

    nperiods = 2
    period = Month(1)
    startdate = DateTime(2016)
    enddate = startdate + period * nperiods

    stack2 = GeoStack(NamedTuple{(:lower, :upper, :temp)}((lowerdata, upperdata, tempdata)))
    stack3 = GeoStack(NamedTuple{(:lower, :upper, :temp)}((lowerdata, upperdata, tempdata)))
    stack4 = GeoStack(NamedTuple{(:lower, :upper, :temp)}((lowerdata, upperdata, tempdata)))

    dimz = (Time([DateTime(2016, 1, 1, 9), DateTime(2016, 1, 16, 15),
                  DateTime(2016, 2, 1, 10), DateTime(2016, 2, 16, 14)]),)
    series = GeoSeries([stack, stack2, stack3, stack4], dimz)
    model = growth, lower, upper

    # TODO test output
    output = mapgrowth(model, series;
                       startdate=startdate,
                       nperiods=4,
                       period=Day(15),
                       subperiod=Day(1),
                       subperiod_starts=subperiod_starts,
                       constructor=identity)
    @test typeof(dims(output)) <: Tuple{Lat,Lon,Time}
    @test length(val(dims(output, Time))) == 4
    dims(output, Time)

    @test_broken output[Time(DateTime(2016, 1, 1))] # == TODO
    @test_broken output[Time(DateTime(2016, 1, 16))] # == TODO
    @test_broken output[Time(DateTime(2016, 1, 31))] # == TODO
    @test_broken output[Time(DateTime(2016, 2, 15))] # == TODO

    output = mapgrowth(model, series;
                       startdate=startdate,
                       nperiods=2,
                       period=Month(1),
                       subperiod=Day(2),
                       subperiod_starts=subperiod_starts,
                       constructor=identity)
    @test length(val(dims(output, Time))) == 2
    @test_broken output[Time(DateTime(2016, 1))] # == TODO
    @test_broken output[Time(DateTime(2016, 2))] # == TODO
end

# We use wget and unzip to handle files, so skip windows
if !Sys.iswindows() 
    @testset "SMAP" begin
        include("smap.jl")
    end
end
