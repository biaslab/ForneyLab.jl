import Base: prod!

export Gaussian, prod!

abstract type Gaussian <: SoftFactor end

function prod!{F1<:Gaussian, F2<:Gaussian}(
    x::ProbabilityDistribution{Univariate, F1},
    y::ProbabilityDistribution{Univariate, F2},
    z::ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}=ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=1.0))

    z.params[:xi] = unsafeWeightedMean(x) + unsafeWeightedMean(y)
    z.params[:w] = unsafePrecision(x) + unsafePrecision(y)

    return z
end

@symmetrical function prod!{F<:Gaussian}(
    x::ProbabilityDistribution{Univariate, F},
    y::ProbabilityDistribution{Univariate, PointMass},
    z::ProbabilityDistribution{Univariate, PointMass}=ProbabilityDistribution(Univariate, PointMass, m=0.0))

    z.params[:m] = y.params[:m]
    return z
end

function prod!{F1<:Gaussian, F2<:Gaussian}(
    x::ProbabilityDistribution{Multivariate, F1},
    y::ProbabilityDistribution{Multivariate, F2},
    z::ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}=ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[NaN], w=[NaN].'))

    z.params[:xi] = unsafeWeightedMean(x) + unsafeWeightedMean(y)
    z.params[:w] = unsafePrecision(x) + unsafePrecision(y)

    return z
end

@symmetrical function prod!{F<:Gaussian}(
    x::ProbabilityDistribution{Multivariate, F},
    y::ProbabilityDistribution{Multivariate, PointMass},
    z::ProbabilityDistribution{Multivariate, PointMass}=ProbabilityDistribution(Multivariate, PointMass, m=[NaN]))

    z.params[:m] = deepcopy(y.params[:m])

    return z
end

function sample{F<:Gaussian}(dist::ProbabilityDistribution{Univariate, F})
    isProper(dist) || error("Cannot sample from improper distribution")
    return sqrt(unsafeCov(dist))*randn() + dist.params[:m]
end

function sample{F<:Gaussian}(dist::ProbabilityDistribution{Multivariate, F})
    isProper(dist) || error("Cannot sample from improper distribution")
    return chol(unsafeCov(dist))' *randn(dims(dist)) + dist.params[:m]
end

# Entropy functional
function differentialEntropy{F<:Gaussian}(dist::ProbabilityDistribution{Univariate, F})
    return  0.5*log(unsafeCov(dist)) +
            0.5*log(2*pi) +
            0.5
end

function differentialEntropy{F<:Gaussian}(dist::ProbabilityDistribution{Multivariate, F})
    return  0.5*log(det(unsafeCov(dist))) +
            (dims(dist)/2)*log(2*pi) +
            (dims(dist)/2)
end