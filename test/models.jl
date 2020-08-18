using GrowthMaps, GeoData, Unitful, Test
using Unitful: °C, K, hr, d, mol, cal
using GrowthMaps: rate, condition, conditionalrate, combinelayers

dimz = Lat((10, 20)), Lon((100, 130))

# Models
lowerdata = [1. 2. 3.
             4. 5. 6.]
upperdata = [9. 8. 7.
             6. 5. 4.]
tempdata = [270.0 280.0 290.0
            300.0 310.0 320.0]

lowerunitful = lowerdata * K
upperunitful = upperdata * K
tempunitful = tempdata * K

threshold = 5K
mortalityrate = -1/K
lower = Layer{:lower,K}(LowerStress(threshold, mortalityrate))

@testset "Lower stress" begin
    @test rate(lower, 1.0K) == -4.
    @test condition(lower, 1.0K) == true
    @test rate(lower, 6.0K) == 1.
    @test condition(lower, 6.0K) == false
    @test keys(lower) == :lower
    @test condition.(Ref(lower), lowerdata * K) == [true true  true
                                                    true false false]
    @test conditionalrate.(Ref(lower), lowerdata) == [-4. -3. -2.
                                                      -1.  0.  0.]
end

threshold = 6K
mortalityrate = -2/K
upper = Layer(:upper, K, UpperStress(threshold, mortalityrate))

@testset "Upper stress" begin
    @test rate(upper, 9.0K) == -6.
    @test condition(upper, 9.0K) == true
    @test rate(upper, 0.0K) == 12.
    @test condition(upper, 0.0K) == false

    @test keys(upper) == :upper
    @test condition.(Ref(upper), upperdata * K) == [true  true  true
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
growth = Layer(:temp, K, SchoolfieldIntrinsicGrowth(p, ΔH_A, ΔH_L, T_halfL, ΔH_H, T_halfH, T_ref))

@testset "Schoolfield Intrinsic growth" begin
    @test keys(growth) == :temp
    @test condition.(Ref(growth), tempdata * K) == Bool[1 1 1
                                                        1 1 1]
    rates = Array(conditionalrate.(Ref(growth), tempdata))
    refrates = [8.034064e-03 3.154654e-02 1.128590e-01; 1.856456e-01 1.007895e-07 7.046076e-14]
    # The original R code was not very accurate, we only use two significant figures.
    @test rates ≈ refrates rtol = 1e-2
end


@testset "Combined layers" begin
    stack = NamedTuple{(:lower, :upper, :temp)}((lowerdata, upperdata, tempdata))
    @test combinelayers((lower, upper), stack) == [-10. -7. -4.
                                                   -1.   0.  0.]
    combinelayers((growth, lower, upper), stack)
end

