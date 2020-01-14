module GrowthMaps

# Use the README as the module docs
@doc read(joinpath(dirname(@__DIR__), "README.md"), String) GrowthMaps

using ConstructionBase,
      Dates,
      FieldMetadata,
      Flatten,
      GeoData,
      LsqFit,
      Interact,
      Plots,
      Unitful,
      UnitfulRecipes

import FieldMetadata: @flattenable, flattenable, @bounds, bounds
using GeoData: Time, rebuild
using Unitful: °C, K, R
using Base: tail

export mapgrowth, fit, manualfit!, mapfit!

export RateModel

export GrowthModel, SchoolfieldIntrinsicGrowth

export StressModel, LowerStress, ColdStress, WiltStress,
       UpperStress, HeatStress, MoistureStress

export Layer, ModelWrapper

include("models.jl")
include("framework.jl")
include("fit.jl")

end # module
