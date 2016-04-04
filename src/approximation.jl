export ApproximationType, Laplace, MomentMatching
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
    Approximation{dist<:ProbabilityDistribution, approx_type<:ApproximationType}

Used to specify an approximation of a probability distribution.
`dist` is the type of the approximating distribution.
`approx_type` is the approximation type. Example:

    Approximation{GaussianDistribution, MomentMatching}
"""
abstract Approximation{dist<:ProbabilityDistribution, approx_type<:ApproximationType}