module GrowthMaps

# Use the README as the module docs
@doc read(joinpath(dirname(@__DIR__), "README.md"), String) GrowthMaps

using ConstructionBase,
      Dates,
      GeoData,
      InteractModels,
      LsqFit,
      Plots,
      Reexport,
      Unitful,
      UnitfulRecipes

@reexport using ModelParameters

using GeoData: rebuild
using Unitful: °C, K, R, Units
using Base: tail

export mapgrowth, fit, fit!, manualfit!, mapfit!

export RateModel

export GrowthModel, SchoolfieldIntrinsicGrowth

export StressModel, LowerStress, UpperStress

export Layer

include("models.jl")
include("framework.jl")
include("fit.jl")

function run_precompile()
    p = Param(3e-1; bounds=(3e-2, 3e0))
    ΔH_A = Param(3e4; units=u"cal/mol", bounds=(3e3, 3e5))
    ΔH_L = Param(-1e5; units=u"cal/mol", bounds=(-1e6, -1e4))
    ΔH_H = Param(3e5; units=u"cal/mol", bounds=(3e4, 3e6))
    Thalf_L = Param(2e2; units=u"K", bounds=(2e1, 2e3))
    Thalf_H = Param(3e2; units=u"K", bounds=(3e1, 3e3))
    T_ref = u"K"(25.0u"°C")
    growthmodel = SchoolfieldIntrinsicGrowth(p, ΔH_A, ΔH_L, Thalf_L, ΔH_H, Thalf_H, T_ref)
    growth = Layer(:x, u"K", growthmodel)
    model = Model(growth)

    coldthresh = Param(ustrip(u"K", -10.0f0u"°C"); units=u"K", bounds=(240, 290))
    coldmort = Param(-log(1.23f0); units=u"K"^-1, bounds=(0, 0.4))
    coldstress = Layer(:x, u"K", LowerStress(coldthresh, coldmort))
    log(1.23f0)

    heatthresh = Param(ustrip(u"K"(30.0u"°C")); units=u"K", bounds=(280, 330))
    heatmort = Param(-log(1.15); units= u"K"^-1, bounds=(0, 0.4))
    heatstress = Layer(:x, u"K", UpperStress(heatthresh, heatmort))
    model = Model((growth, coldstress, heatstress))

    return nothing
end

run_precompile()

end # module
