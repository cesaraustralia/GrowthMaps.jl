using GrowthMaps, GeoData, HDF5, Dates, Unitful, Test
using GeoData: rebuild
using Unitful: °C, K, hr, d, mol, cal
using GrowthMaps: rate, condition, conditionalrate, combinemodels, periodstartdates

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
        @test periodstartdates(startdate, period, nperiods) == [DateTime(2016, 1, 1), DateTime(2016, 2, 1)]
    end
end

@testset "Integration" begin
    dimz = Lat(10:10), Lon((100, 110))

    # Set up series data
    stressdata = GeoArray.([[1. 2.], [1. 2.],
                            [2. 3.], [3. 4.], [2.5 3.5],
                            [4. 5.], [5. 6.],
                            [6. 7.], [6. 7.], [6. 7.]], Ref(dimz); name="stress")
    tempdata = GeoArray.([[270. 280.], [270. 280.],
                          [270. 280.], [270. 280.], [270. 280.],
                          [270. 280.], [270. 280.],
                          [270. 280.], [270. 280.], [270. 280.]], Ref(dimz); name="tempdata")

    # Build a GeoSeries
    stacks = [GeoStack(NamedTuple{(:stress, :temp)}((stressdata[i], tempdata[i]))) for i in 1:length(stressdata)]
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
    @test series[At(DateTime(2016, 1, 3, 9))][:stress] == [1. 2.]

    # Set up models
    lowerthreshold = 5K; lowermortalityrate = -1/K
    upperthreshold = 5K; uppermortalityrate = -1/K
    lower = Layer(:stress, LowerStress(lowerthreshold, lowermortalityrate))
    upper = Layer(:stress, UpperStress(upperthreshold, uppermortalityrate))

    # Lower
    output = mapgrowth(lower, series;
        period=Month(1),
        nperiods=4,
        startdate=DateTime(2016, 1, 3),
    );

    # Test were are not touching the original arrays
    @test series[At(DateTime(2016, 1, 3, 9))][:stress] == [1. 2.]

    @test output[Ti(1)] == [-4.0 -3.0]
    @test output[Ti(2)] == [-2.5 -1.5]
    @test output[Ti(At(DateTime(2016, 3, 3)))] == [-0.5 0.0]
    @test output[Ti(At(DateTime(2016, 4, 3)))] == [0.0 0.0]

    @test typeof(dims(output)) <: Tuple{Lat,Lon,Ti}
    @test length(val(dims(output, Ti))) == 4

    # Upper
    output = mapgrowth(upper, series;
        period=Month(1),
        nperiods=4,
        startdate=DateTime(2016, 1, 3),
    );

    @test output[Ti(1)] == [0. 0.]
    @test output[Ti(2)] == [0. 0.]
    @test output[Ti(At(DateTime(2016, 3, 3)))] == [0.0 -0.5 ]
    @test output[Ti(At(DateTime(2016, 4, 3)))] == [-1.0 -2.0]

    # Lower and Upper
    output = mapgrowth((lower, upper), series;
        period=Month(1),
        nperiods=4,
        startdate=DateTime(2016, 1, 3),
    );

    # Test were are still not touching the original arrays
    @test series[At(DateTime(2016, 1, 3, 9))][:stress] == [1. 2.]

    @test output[Ti(1)] == [-4.0 -3.0]
    @test output[Ti(2)] == [-2.5 -1.5]
    @test output[Ti(At(DateTime(2016, 3, 3)))] == [-0.5 -0.5]
    @test output[Ti(At(DateTime(2016, 4, 3)))] == [-1.0 -2.0]

    # TODO test temp
end
