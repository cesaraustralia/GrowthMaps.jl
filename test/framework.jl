using GrowthMaps, Rasters, ModelParameters, Dates, Test, StaticArrays
using Unitful: °C, K, hr, d, mol, cal

dimz = Y(10:10), X(100:10:110)

# Set up series data
stressdata = Raster.([[1. 2.], [1. 2.],
                      [2. 3.], [3. 4.], [2.5 3.5],
                      [4. 5.], [5. 6.],
                      [6. 7.], [6. 7.], [6. 7.]], Ref(dimz); name=:stress)
# TODO set up tempdata
tempdata = Raster.([[270. 280.], [270. 280.],
                    [270. 280.], [270. 280.], [270. 280.],
                    [270. 280.], [270. 280.],
                    [270. 280.], [270. 280.], [270. 280.]], Ref(dimz); name=:tempdata)

# Build a RasterSeries
stacks = [RasterStack(NamedTuple{(:stress, :temp)}((stressdata[i], tempdata[i]))) for i in 1:length(stressdata)]
tspan = [
    DateTime(2016, 1, 3, 9),
    DateTime(2016, 1, 6, 15),
    DateTime(2016, 2, 3, 10),
    DateTime(2016, 2, 3, 14),
    DateTime(2016, 2, 18, 10),
    DateTime(2016, 3, 3, 3),
    DateTime(2016, 3, 3, 8),
    DateTime(2016, 4, 3, 14),
    DateTime(2016, 4, 4, 10),
    DateTime(2016, 4, 16, 14)
]; 
timedim = Ti(tspan; span=Rasters.Regular(Hour(3)))
series = RasterSeries(stacks, timedim)

@test series[At(DateTime(2016, 1, 3, 9))][:stress] == [1. 2.]

@testset "mapgrowth float" begin
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

    @test typeof(dims(output)) <: Tuple{Y,X,Ti}
    @test length(dims(output, Ti)) == 4

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

struct MultiRate <: RateModel end

@inline GrowthMaps.rate(l::MultiRate, x) = SA[x/K, 2x/K, 3x/K]

@testset "mapgrowth returning static array" begin
    model = Layer(:stress, K, MultiRate())

    output = mapgrowth(model; 
        series=series,
        tspan=DateTime(2016, 1, 3):Month(1):DateTime(2016, 4, 3),
        initval=SA[0.0, 0.0, 0.0]
    );
    @test output[1, 1, 1] == SA[1.0, 2.0, 3.0]
    @test output[1, 2, 1] == SA[2.0, 4.0, 6.0]
end

