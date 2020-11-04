using GrowthMaps, GeoData, ModelParameters, Dates, Test
using Unitful: °C, K, hr, d, mol, cal

@testset "mapgrowth" begin
    dimz = Lat(10:10), Lon((100, 110))

    # Set up series data
    stressdata = GeoArray.([[1. 2.], [1. 2.],
                            [2. 3.], [3. 4.], [2.5 3.5],
                            [4. 5.], [5. 6.],
                            [6. 7.], [6. 7.], [6. 7.]], Ref(dimz); name=:stress)
    # TODO set up tempdata
    tempdata = GeoArray.([[270. 280.], [270. 280.],
                          [270. 280.], [270. 280.], [270. 280.],
                          [270. 280.], [270. 280.],
                          [270. 280.], [270. 280.], [270. 280.]], Ref(dimz); name=:tempdata)

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
    lowerthreshold = 5
    lowermortalityrate = -1
    lower = Layer(:stress, LowerStress(lowerthreshold, lowermortalityrate))

    upperthreshold = 5K
    uppermortalityrate = -1/K
    upper = Layer(:stress, K, UpperStress(upperthreshold, uppermortalityrate))

    # Lower
    output = mapgrowth(lower; 
        series=series,
        tspan=DateTime(2016, 1, 3):Month(1):DateTime(2016, 4, 3)
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
    output = mapgrowth(upper;
        series=series,
        tspan=DateTime(2016, 1, 3):Month(1):DateTime(2016, 4, 3)
    );

    @test output[Ti(1)] == [0. 0.]
    @test output[Ti(2)] == [0. 0.]
    @test output[Ti(At(DateTime(2016, 3, 3)))] == [0.0 -0.5 ]
    @test output[Ti(At(DateTime(2016, 4, 3)))] == [-1.0 -2.0]

    # Lower and Uppera, in a Model with an extra period
    output = mapgrowth(Model((lower, upper));
        series=series,
        tspan=DateTime(2016, 1, 3):Month(1):DateTime(2016, 5, 3)
    );

    # Test were are still not touching the original arrays
    @test series[At(DateTime(2016, 1, 3, 9))][:stress] == [1. 2.]

    @test output[Ti(1)] == [-4.0 -3.0]
    @test output[Ti(2)] == [-2.5 -1.5]
    @test output[Ti(At(DateTime(2016, 3, 3)))] == [-0.5 -0.5]
    @test output[Ti(At(DateTime(2016, 4, 3)))] == [-1.0 -2.0]

    @test_logs (:warn,"No files found for the 1 month period starting 2016-05-03T00:00:00") 
    mapgrowth(
        Model((lower, upper));
        series=series,
        tspan=DateTime(2016, 1, 3):Month(1):DateTime(2016, 5, 3)
    );

end
