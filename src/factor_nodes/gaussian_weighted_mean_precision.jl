export GaussianWeightedMeanPrecision

"""
Description:

    A Gaussian with weighted-mean-precision parameterization:

    f(x,xi,w) = 𝒩(x | xi/w, 1/w)

Construction:

    No node constructor available
"""
abstract GaussianWeightedMeanPrecision <: Gaussian

ProbabilityDistribution(::Type{GaussianWeightedMeanPrecision}) = ProbabilityDistribution(GaussianWeightedMeanPrecision, xi=0.0, w=1.0)

function prod!{T<:Gaussian, U<:Gaussian}(   x::ProbabilityDistribution{T},
                                            y::ProbabilityDistribution{U},
                                            z::ProbabilityDistribution{GaussianWeightedMeanPrecision}=ProbabilityDistribution(GaussianWeightedMeanPrecision))

    x = convert(ProbabilityDistribution{GaussianWeightedMeanPrecision}, x)
    y = convert(ProbabilityDistribution{GaussianWeightedMeanPrecision}, y)

    z.params[:xi] = x.params[:xi] + y.params[:xi]
    z.params[:w] = x.params[:w] + y.params[:w]

    return z
end

convert(::Type{ProbabilityDistribution{GaussianWeightedMeanPrecision}}, dist::ProbabilityDistribution{GaussianMeanVariance}) =
    ProbabilityDistribution(GaussianWeightedMeanPrecision, xi=dist.params[:m]/dist.params[:v], w=1/dist.params[:v])

convert(::Type{ProbabilityDistribution{GaussianMeanVariance}}, dist::ProbabilityDistribution{GaussianWeightedMeanPrecision}) =
    ProbabilityDistribution(GaussianMeanVariance, m=dist.params[:xi]/dist.params[:w], v=1/dist.params[:w])