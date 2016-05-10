export ApproximationType, Laplace, MomentMatching, ExponentMomentMatching, LogMomentMatching
export Approximation

"""
`ApproximationType` has subtypes that represent approximation types, such as `Laplace` and `MomentMatching`.
"""
abstract ApproximationType

"""
`Laplace <: ApproximationType` represents the Gaussian approximation through Laplace's method.
"""
abstract Laplace <: ApproximationType

"""
`MomentMatching <: ApproximationType` represents the moment matching approximation.
"""
abstract MomentMatching <: ApproximationType

"""
`ExpMomentMatching <: ApproximationType` represents the approximation that matches the moments of `exp(X)`, where X is the random variable.
"""
abstract ExponentMomentMatching <: ApproximationType

"""
`LogMomentMatching <: ApproximationType` represents the approximation that matches the moments of `ln(X)`, where X is the random variable.
"""
abstract LogMomentMatching <: ApproximationType

"""
    Approximation{dist<:ProbabilityDistribution, approx_type<:ApproximationType}

Used to specify an approximation of a probability distribution.
`dist` is the type of the approximating distribution.
`approx_type` is the approximation type. Example:

    Approximation{Gaussian, MomentMatching}
"""
abstract Approximation{dist<:ProbabilityDistribution, approx_type<:ApproximationType}