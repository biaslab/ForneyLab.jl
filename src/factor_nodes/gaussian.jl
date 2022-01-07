export Gaussian, prod!, convert

# Convert parameterizations
function convert(::Type{ProbabilityDistribution{V, GaussianMeanPrecision}}, dist::ProbabilityDistribution{V, GaussianMeanVariance}) where V<:VariateType
    w = cholinv(dist.params[:v])
    m = deepcopy(dist.params[:m])

    return ProbabilityDistribution(V, GaussianMeanPrecision, m=m, w=w)
end

function convert(::Type{ProbabilityDistribution{V, GaussianMeanPrecision}}, dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType
    w = deepcopy(dist.params[:w])
    m = cholinv(w)*dist.params[:xi]

    return ProbabilityDistribution(V, GaussianMeanPrecision, m=m, w=w)
end

function convert(::Type{ProbabilityDistribution{V, GaussianMeanVariance}}, dist::ProbabilityDistribution{V, GaussianMeanPrecision}) where V<:VariateType
    v = cholinv(dist.params[:w])
    m = deepcopy(dist.params[:m])

    return ProbabilityDistribution(V, GaussianMeanVariance, m=m, v=v)
end

function convert(::Type{ProbabilityDistribution{V, GaussianMeanVariance}}, dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType
    v = cholinv(dist.params[:w])
    m = v*dist.params[:xi]

    return ProbabilityDistribution(V, GaussianMeanVariance, m=m, v=v)
end

function convert(::Type{ProbabilityDistribution{V, GaussianWeightedMeanPrecision}}, dist::ProbabilityDistribution{V, GaussianMeanPrecision}) where V<:VariateType
    w = deepcopy(dist.params[:w])
    xi = w*dist.params[:m]

    return ProbabilityDistribution(V, GaussianWeightedMeanPrecision, xi=xi, w=w)
end

function convert(::Type{ProbabilityDistribution{V, GaussianWeightedMeanPrecision}}, dist::ProbabilityDistribution{V, GaussianMeanVariance}) where V<:VariateType
    w = cholinv(dist.params[:v])
    xi = w*dist.params[:m]

    return ProbabilityDistribution(V, GaussianWeightedMeanPrecision, xi=xi, w=w)
end

# Convert VariateTypes
convert(::Type{ProbabilityDistribution{Multivariate, GaussianMeanVariance}}, dist::ProbabilityDistribution{Univariate, GaussianMeanVariance}) =
    ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[dist.params[:m]], v=mat(dist.params[:v]))
convert(::Type{ProbabilityDistribution{Multivariate, GaussianMeanPrecision}}, dist::ProbabilityDistribution{Univariate, GaussianMeanPrecision}) =
    ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=[dist.params[:m]], w=mat(dist.params[:w]))
convert(::Type{ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}}, dist::ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}) =
    ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[dist.params[:xi]], w=mat(dist.params[:w]))

function prod!(
    x::ProbabilityDistribution{Univariate, <:Gaussian},
    y::ProbabilityDistribution{Univariate, <:Gaussian},
    z::ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}=ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=1.0))

    z.params[:xi] = unsafeWeightedMean(x) + unsafeWeightedMean(y)
    z.params[:w] = unsafePrecision(x) + unsafePrecision(y)

    return z
end

@symmetrical function prod!(
    x::ProbabilityDistribution{Univariate, <:Gaussian},
    y::ProbabilityDistribution{Univariate, PointMass},
    z::ProbabilityDistribution{Univariate, PointMass}=ProbabilityDistribution(Univariate, PointMass, m=0.0))

    z.params[:m] = y.params[:m]
    return z
end

function prod!(
    x::ProbabilityDistribution{Multivariate, <:Gaussian},
    y::ProbabilityDistribution{Multivariate, <:Gaussian},
    z::ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}=ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[NaN], w=transpose([NaN])))

    z.params[:xi] = unsafeWeightedMean(x) + unsafeWeightedMean(y)
    z.params[:w] = unsafePrecision(x) + unsafePrecision(y)

    return z
end

@symmetrical function prod!(
    x::ProbabilityDistribution{Multivariate, <:Gaussian},
    y::ProbabilityDistribution{Multivariate, PointMass},
    z::ProbabilityDistribution{Multivariate, PointMass}=ProbabilityDistribution(Multivariate, PointMass, m=[NaN]))

    z.params[:m] = deepcopy(y.params[:m])

    return z
end

function sample(dist::ProbabilityDistribution{Univariate, <:Gaussian})
    isProper(dist) || error("Cannot sample from improper distribution")
    (m,v) = unsafeMeanCov(dist)
    return sqrt(v)*randn() + m
end

function sample(dist::ProbabilityDistribution{Univariate, <:Gaussian}, n_samples::Int64)
    isProper(dist) || error("Cannot sample from improper distribution")
    (m,v) = unsafeMeanCov(dist)

    return sqrt(v).*randn(n_samples) .+ m
end

function sample(dist::ProbabilityDistribution{Multivariate, <:Gaussian})
    isProper(dist) || error("Cannot sample from improper distribution")
    (m,V) = unsafeMeanCov(dist)
    return (cholesky(default_cholesky_mode, V)).U' *randn(dims(dist)) + m
end

function sample(dist::ProbabilityDistribution{Multivariate, <:Gaussian}, n_samples::Int64)
    isProper(dist) || error("Cannot sample from improper distribution")
    (m,V) = unsafeMeanCov(dist)
    U = (cholesky(default_cholesky_mode, V)).U
    d = dims(dist)

    return [U' *randn(d) + m for i in 1:n_samples]
end

function naturalParams(dist::ProbabilityDistribution{<:VariateType, <:Gaussian})
    (xi, w) = unsafeWeightedMeanPrecision(dist)
    return [xi, -0.5*w] # Returns vector and matrix for Multivariate to retain dimensionality info; η[2] is not positive definite
end

standardDistribution(::Type{V}, ::Type{<:Gaussian}; η::Vector) where V<:VariateType = ProbabilityDistribution(V, GaussianWeightedMeanPrecision, xi=η[1], w=-2*η[2])

logNormalizer(::Type{Univariate}, ::Type{<:Gaussian}; η::Vector) = -η[1]^2/(4*η[2]) - 0.5*log(-2*η[2])
logNormalizer(::Type{Multivariate}, ::Type{<:Gaussian}; η::Vector) = η[1]'*cholinv(-4*η[2])*η[1] - 0.5*logdet(-2*η[2])

logPdf(V::Type{Univariate}, ::Type{F}, x::Number; η::Vector) where F<:Gaussian = -0.5*log(2pi) + [x, x^2]'*η - logNormalizer(V, F; η=η)
logPdf(V::Type{Multivariate}, ::Type{F}, x::Vector; η::Vector) where F<:Gaussian = -0.5*length(η[1])*log(2pi) + [x; vec(x*x')]'*[η[1]; vec(η[2])] - logNormalizer(V, F; η=η)

# Entropy functional
function differentialEntropy(dist::ProbabilityDistribution{Univariate, <:Gaussian})
    return  0.5*log(unsafeCov(dist)) +
            0.5*log(2*pi) +
            0.5
end

function differentialEntropy(dist::ProbabilityDistribution{Multivariate, <:Gaussian})
    d = dims(dist)[1]
    return  0.5*logdet(unsafeCov(dist)) +
            (d/2)*log(2*pi) +
            d/2
end
