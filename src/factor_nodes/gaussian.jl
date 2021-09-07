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
    z::ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}=ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[NaN], w=transpose([NaN]))) where {F1<:Gaussian, F2<:Gaussian}

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

function sample(dist::ProbabilityDistribution{Univariate, F}, n_samples::Int64) where F<:Gaussian
    isProper(dist) || error("Cannot sample from improper distribution")
    (m,v) = unsafeMeanCov(dist)

    return sqrt(v).*randn(n_samples) .+ m
end

function sample(dist::ProbabilityDistribution{Multivariate, F}) where F<:Gaussian
    isProper(dist) || error("Cannot sample from improper distribution")
    (m,V) = unsafeMeanCov(dist)
    return (cholesky(Matrix(Hermitian(V)))).U' *randn(dims(dist)) + m
end

function sample(dist::ProbabilityDistribution{Multivariate, F}, n_samples::Int64) where F<:Gaussian
    isProper(dist) || error("Cannot sample from improper distribution")
    (m,V) = unsafeMeanCov(dist)
    U = (cholesky(Matrix(Hermitian(V)))).U
    d = dims(dist)

    return [U' *randn(d) + m for i in 1:n_samples]
end

# Standard parameters to natural parameters
naturalParams(dist::ProbabilityDistribution{Univariate, F}) where F<:Gaussian = [unsafeWeightedMean(dist), -0.5*unsafePrecision(dist)]

naturalParams(dist::ProbabilityDistribution{Multivariate, F}) where F<:Gaussian = [vec(unsafeWeightedMean(dist)); vec(-0.5*unsafePrecision(dist))]

# Natural parameters to standard dist. type
function standardDist(dist::ProbabilityDistribution{Univariate, F}, η::Vector) where F<:Gaussian
    ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision,xi=η[1],w=-2*η[2])
end

function standardDist(dist::ProbabilityDistribution{Multivariate, F}, η::Vector) where F<:Gaussian
    d = dims(dist)
    XI, W = η[1:d], reshape(-2*η[d+1:end],d,d)
    W = Matrix(Hermitian(W + tiny*diageye(d))) # Ensure precision is always invertible
    ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision,xi=XI,w=W)
end

# Natural parameters to standard message type
function standardMessage(dist::ProbabilityDistribution{Univariate, F}, η::Vector) where F<:Gaussian
    Message(Univariate, GaussianWeightedMeanPrecision,xi=η[1],w=-2*η[2])
end

function standardMessage(dist::ProbabilityDistribution{Multivariate, F}, η::Vector) where F<:Gaussian
    d = dims(dist)
    XI, W = η[1:d], reshape(-2*η[d+1:end],d,d)
    W = Matrix(Hermitian(W + tiny*diageye(d))) # Ensure precision is always invertible
    Message(Multivariate, GaussianWeightedMeanPrecision,xi=XI,w=W)
end

function logNormalizer(dist::ProbabilityDistribution{Univariate, F}, η::Vector) where F<:Gaussian
    return - η[1]^2/(4*η[2]) - 0.5*log(-2*η[2])
end

function logNormalizer(dist::ProbabilityDistribution{Multivariate, F}, η::Vector) where F<:Gaussian
    d = dims(dist)
    return - 0.25*η[1:d]'*(reshape(η[d+1:end],(d,d))\η[1:d]) -0.5*logdet(-2*reshape(η[d+1:end],(d,d)))
end

# logPdf wrt natural params. ForwardDiff is not stable with reshape function which
# precludes the usage of logPdf functions previously defined. Below function is
# meant to be used with Zygote.
function logPdf(dist::ProbabilityDistribution{Univariate, F}, η::Vector, x) where F<:Gaussian
    h(x) = 1/sqrt(2*pi)
    ϕ(x) = [x,x^2]
    return log(h(x)) + transpose(ϕ(x))*η - logNormalizer(dist,η)
end

function logPdf(dist::ProbabilityDistribution{Multivariate, F}, η::Vector, x) where F<:Gaussian
    d = dims(dist)
    h(x) = (2*pi)^(-0.5*d)
    ϕ(x) = [x;vec(x*transpose(x))]
    return log(h(x)) + transpose(ϕ(x))*η - logNormalizer(dist,η)
end

# Entropy functional
function differentialEntropy(dist::ProbabilityDistribution{Univariate, F}) where F<:Gaussian
    return  0.5*log(unsafeCov(dist)) +
            0.5*log(2*pi) +
            0.5
end

function differentialEntropy(dist::ProbabilityDistribution{Multivariate, F}) where F<:Gaussian
    return  0.5*log(det(unsafeCov(dist))) +
            (dims(dist)/2)*log(2*pi) +
            (dims(dist)/2)
end
