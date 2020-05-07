using GrowthMaps, GeoData, Unitful, Test, Plots
using Unitful: °C, K, hr, d, mol, cal
using GrowthMaps: rate, condition, conditionalrate, combinemodels

#= This file just tests that the fitting routines run, not
that they work, which is difficult with Interact.jl. =#

p = 0.3
ΔH_A = 2e4cal * mol^-1
ΔH_L = -1e5cal * mol^-1
ΔH_H = 3e5cal * mol^-1
T_halfL = 250.0K
T_halfH = 300.0K
T_ref = K(25.0°C)
growth = Layer(:temp, K, SchoolfieldIntrinsicGrowth(p, ΔH_A, ΔH_L, T_halfL, ΔH_H, T_halfH, T_ref))
wrapper = ModelWrapper(growth)

obs = [(10.0, 0.0), (20.0, 0.2), (30.0, 0.1), (40.0, 0.0)]

@testset "Auto fit" begin
    fitted = fit!(wrapper, obs)
end

@testset "Manual fit interface" begin
    range = NamedTuple{(:temp,)}((0°C:1°C:45°C,))
    interface = manualfit!(wrapper, range; obs=obs);
end


@testset "Map fit interface" begin
    ## Set up series data
    tempdata = GeoArray.([rand(0.0:40.0, 10, 10) for i in 1:5], Ref((Lon, Lat)); 
                         name="stress", missingval=-99.0)

    stacks = [GeoStack(NamedTuple{(:temp,)}((tempdata[i],))) for i in 1:length(tempdata)]
    timedim = (Ti((1:1:5)hr; mode=Sampled(; span=Regular(1hr))),)
    modelkwargs = (series=GeoSeries(stacks, timedim), period=2hr, nperiods=3, startdate=1hr)
    
    mapfit!(wrapper, modelkwargs; occurrence=[(1, 2), (9, 10)]);
end

nothing
