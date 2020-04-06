using GrowthMaps, GeoData, HDF5, Dates, Unitful, Test
using GeoData: Time, rebuild
using Unitful: °C, K, hr, d, mol, cal
using GrowthMaps: rate, condition, conditionalrate, combinemodels, 
      period_startdates, subset_startdates

dimz = Lat((10, 20)), Lon((100, 130))

# Models
lowerdata = GeoArray([1. 2. 3.
                      4. 5. 6.], dimz)
upperdata = GeoArray([9. 8. 7.
                      6. 5. 4.], dimz)
tempdata = GeoArray([270.0 280.0 290.0
                     300.0 310.0 320.0], dimz)

lowerunitful = lowerdata * K
upperunitful = upperdata * K
tempunitful = tempdata * K

threshold = 5K
mortalityrate = -1/K
lower = Layer(:lower, LowerStress(threshold, mortalityrate))

@testset "Lower stress" begin
    @test keys(lower) == :lower
    @test condition.(Ref(lower), lowerdata) == [true true  true
                                                true false false]
    @test conditionalrate.(Ref(lower), lowerdata) == [-4. -3. -2.
                                                      -1.  0.  0.]
end

threshold = 6K
mortalityrate = -2/K
upper = Layer(:upper, UpperStress(threshold, mortalityrate))

@testset "Upper stress" begin
    @test keys(upper) == :upper
    @test condition.(Ref(upper), upperdata) == [true  true  true
                                                false false false]
    @test conditionalrate.(Ref(upper), upperdata) == [-6. -4. -2.
                                                       0.  0.  0.]
end

p = 0.3
ΔH_A = 2e4cal * mol^-1
ΔH_L = -1e5cal * mol^-1
ΔH_H = 3e5cal * mol^-1
T_halfL = 250.0K
T_halfH = 300.0K
T_ref = K(25.0°C)
growth = Layer(:temp, SchoolfieldIntrinsicGrowth(p, ΔH_A, ΔH_L, T_halfL, ΔH_H, T_halfH, T_ref))

@testset "Schoolfield Intrinsic growth" begin
    @test keys(growth) == :temp
    @test_broken condition.(Ref(growth), tempdata) == [-0. -0. -0.
                                                       -0.  0.  0.]
    rates = Array(conditionalrate.(Ref(growth), tempdata))
    refrates = [8.034064e-03 3.154654e-02 1.128590e-01; 1.856456e-01 1.007895e-07 7.046076e-14]
    # The original R code was not very accurate, we only use two significant figures.
    @test rates ≈ refrates rtol = 1e-2
end


@testset "Combined models" begin
    stack = GeoStack(NamedTuple{(:lower, :upper, :temp)}((lowerdata, upperdata, tempdata)))
    @test combinemodels((lower, upper), stack) == [-10. -7. -4.
                                                   -1.   0.  0.]
    combinemodels((growth, lower, upper), stack)
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
    nperiods = 4
    period = Month(1)
    subperiod = Day(1)
    startdate = DateTime(2016, 1, 3)
    enddate = startdate + period * nperiods
    subperiod_starts = DateTime.(2016, [1,2,3,4], [3, 3, 3, 3])
    dimz = Lat(10:10), Lon((100, 110))

    lowerdata = GeoArray.([[1. 2.], [-100. -100.], 
                           [2. 3.], [3. 4.], [-100. -100.], 
                           [4. 5.], [5. 6.], 
                           [6. 7.], [-100. -100.], [-100. -100.]], Ref(dimz); name="lower")
    upperdata = GeoArray.([[1. 8.], [9. 8.], 
                           [9. 8.], [9. 8.], [9. 8.], 
                           [9. 8.], [9. 8.], 
                           [9. 8.], [9. 8.], [9. 8.]], Ref(dimz); name="upper")
    tempdata = GeoArray.([[270.0 280.0], [270.0 280.0], 
                          [270.0 280.0], [270.0 280.0], [270.0 280.0],
                          [270.0 280.0], [270.0 280.0], 
                          [270.0 280.0], [270.0 280.0], [270.0 280.0]], Ref(dimz); name="tempdata")
    stacks = [GeoStack(NamedTuple{(:lower, :upper, :temp)}((lowerdata[i], upperdata[i], tempdata[i]))) for i in 1:length(lowerdata)]
    parent(copy(stacks)[5][:lower]) === parent(stacks[5][:lower])
    timedim = (Ti([DateTime(2016, 1, 3, 9), 
                   DateTime(2016, 1, 6, 15),
                   DateTime(2016, 2, 3, 10), 
                   DateTime(2016, 2, 3, 14),
                   DateTime(2016, 2, 18, 10), 
                   DateTime(2016, 3, 3, 3),
                   DateTime(2016, 3, 3, 8), 
                   DateTime(2016, 4, 3, 14),
                   DateTime(2016, 4, 4, 10), 
                   DateTime(2016, 4, 16, 14)
                  ]; mode=Sampled(; span=Regular(Hour(3)))),)

    series = GeoSeries(stacks, timedim)
    threshold = 5K; mortalityrate = -1/K
    lower = Layer(:lower, LowerStress(threshold, mortalityrate))
    model = lower
    output = mapgrowth(model, series;
                       startdate=startdate,
                       nperiods=nperiods,
                       period=period,
                       subperiod=subperiod,
                       subperiod_starts=subperiod_starts)

    @test output[Ti(1)] == [-4.0 -3.0]
    @test output[Ti(2)] == [-2.5 -1.5]
    @test output[Ti(At(DateTime(2016, 3, 3)))] == [-0.5 0.0]
    @test output[Ti(At(DateTime(2016, 4, 3)))] == [0.0 0.0]

    @test typeof(dims(output)) <: Tuple{Lat,Lon,Ti}
    @test length(val(dims(output, Ti))) == 4
end
