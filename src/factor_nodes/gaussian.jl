import Base: prod!, convert

export Gaussian, prod!, convert

function convert(::Type{ProbabilityDistribution{V, GaussianMeanPrecision}},
                 dist::ProbabilityDistribution{V, GaussianMeanVariance}) where V<:VariateType
    w = inv(dist.params[:v])
    m = deepcopy(dist.params[:m])

    return ProbabilityDistribution(V, GaussianMeanPrecision, m=m, w=w)
end

function convert(::Type{ProbabilityDistribution{V, GaussianMeanPrecision}},
                 dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType
    w = deepcopy(dist.params[:w])
    m = inv(w)*dist.params[:xi]

    return ProbabilityDistribution(V, GaussianMeanPrecision, m=m, w=w)
end

function convert(::Type{ProbabilityDistribution{V, GaussianMeanVariance}},
                 dist::ProbabilityDistribution{V, GaussianMeanPrecision}) where V<:VariateType
    v = inv(dist.params[:w])
    m = deepcopy(dist.params[:m])

    return ProbabilityDistribution(V, GaussianMeanVariance, m=m, v=v)
end

function convert(::Type{ProbabilityDistribution{V, GaussianMeanVariance}},
                 dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType
    v = inv(dist.params[:w])
    m = v*dist.params[:xi]

    return ProbabilityDistribution(V, GaussianMeanVariance, m=m, v=v)
end

function convert(::Type{ProbabilityDistribution{V, GaussianWeightedMeanPrecision}},
                 dist::ProbabilityDistribution{V, GaussianMeanPrecision}) where V<:VariateType
    w = deepcopy(dist.params[:w])
    xi = w*dist.params[:m]

    return ProbabilityDistribution(V, GaussianWeightedMeanPrecision, xi=xi, w=w)
end

function convert(::Type{ProbabilityDistribution{V, GaussianWeightedMeanPrecision}},
                 dist::ProbabilityDistribution{V, GaussianMeanVariance}) where V<:VariateType
    w = inv(dist.params[:v])
    xi = w*dist.params[:m]

    return ProbabilityDistribution(V, GaussianWeightedMeanPrecision, xi=xi, w=w)
end

function prod!(
    x::ProbabilityDistribution{Univariate, F1},
    y::ProbabilityDistribution{Univariate, F2},
    z::ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}=ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=1.0)) where {F1<:Gaussian, F2<:Gaussian}

    z.params[:xi] = unsafeWeightedMean(x) + unsafeWeightedMean(y)
    z.params[:w] = unsafePrecision(x) + unsafePrecision(y)

    return z
end

@symmetrical function prod!(
    x::ProbabilityDistribution{Univariate, F},
    y::ProbabilityDistribution{Univariate, PointMass},
    z::ProbabilityDistribution{Univariate, PointMass}=ProbabilityDistribution(Univariate, PointMass, m=0.0)) where F<:Gaussian

    z.params[:m] = y.params[:m]
    return z
end

function prod!(
    x::ProbabilityDistribution{Multivariate, F1},
    y::ProbabilityDistribution{Multivariate, F2},
    z::ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}=ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[0], w=diageye(1))) where {F1<:Gaussian, F2<:Gaussian}

    z.params[:xi] = unsafeWeightedMean(x) + unsafeWeightedMean(y)
    z.params[:w] = unsafePrecision(x) + unsafePrecision(y)

    return z
end

@symmetrical function prod!(
    x::ProbabilityDistribution{Multivariate, F},
    y::ProbabilityDistribution{Multivariate, PointMass},
    z::ProbabilityDistribution{Multivariate, PointMass}=ProbabilityDistribution(Multivariate, PointMass, m=[NaN])) where F<:Gaussian

    z.params[:m] = deepcopy(y.params[:m])

    return z
end

function sample(dist::ProbabilityDistribution{Univariate, F}) where F<:Gaussian
    isProper(dist) || error("Cannot sample from improper distribution")
    (m,v) = unsafeMeanCov(dist)
    return sqrt(v)*randn() + m
end

function sample(dist::ProbabilityDistribution{Multivariate, F}) where F<:Gaussian
    isProper(dist) || error("Cannot sample from improper distribution")
    (m,V) = unsafeMeanCov(dist)
    if isa(V, PDMat)
        return (V.chol).U' *randn(dims(dist)) + m
    else
        return cholesky(V).U' *randn(dims(dist)) + m
    end
end

# Entropy functional
function differentialEntropy(dist::ProbabilityDistribution{Univariate, F}) where F<:Gaussian
    return  0.5*log(unsafeCov(dist)) +
            0.5*log(2*pi) +
            0.5
end

function differentialEntropy(dist::ProbabilityDistribution{Multivariate, GaussianMeanVariance})
    return  0.5*logdet(dist.params[:v]) +
            (dims(dist)/2)*log(2*pi) +
            (dims(dist)/2)
end

function differentialEntropy(dist::ProbabilityDistribution{Multivariate, GaussianMeanPrecision})
    return  0.5*logdet(inv(dist.params[:w])) +
            (dims(dist)/2)*log(2*pi) +
            (dims(dist)/2)
end

function differentialEntropy(dist::ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision})
    return  0.5*logdet(inv(dist.params[:w])) +
            (dims(dist)/2)*log(2*pi) +
            (dims(dist)/2)
end
