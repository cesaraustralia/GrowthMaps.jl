"""
Passed to an optimiser to facilitate fitting any `RateModel`, 
without boilerplate or methods rewrites.
"""
struct ModelFunctor{M<:RateModel}
    model::M
end

"""
Update the model with passed in params, and run [`rate`] over the independent variables.
"""
(f::ModelFunctor)(xs, params) = 
    GrowthMaps.conditionalrate.(Ref(reconstruct(f.model, params, Real)), xs)

"""
    fitlsq(model::RateModel, xs::AbstractArray, ys::AbstractArray)

Fit a model to data with least squares regression, using `curve_fit` from
LsqFit.jl. The passed in model should be initialised with sensible defaults,
these will be used as the initial parameters for the optimization.

Any (nested) `Real` fields on the struct are flattened to a parameter vector using 
[Flatten.jl](http://github.com/rafaqz/Flatten.jl). Fields can be marked to ignore 
using the `@flattenable` macro from [FieldMetadata.jl](http://github.com/rafaqz/FieldMetadata.jl).

## Arguments
- `model`: Any constructed [`RateModel`](@ref), including custom models.
- `xs`: AbstactArray of independent variables (model input values)
- `ys`: AbstactArray of dependent variables (model output values)

## Returns
The model with fitted parameters
"""
fit(model::RateModel, xs::AbstractArray, ys::AbstractArray) = begin
    fit = curve_fit(ModelFunctor(model), xs, ys, [flatten(model, Real)...])
    reconstruct(model, fit.param)
end
