export Gaussian, prod!, convert

# Convert parameterizations
function convert(::Type{Distribution{V, GaussianMeanPrecision}}, dist::Distribution{V, GaussianMeanVariance}) where V<:VariateType
    w = cholinv(dist.params[:v])
    m = deepcopy(dist.params[:m])

    return Distribution(V, GaussianMeanPrecision, m=m, w=w)
end

function convert(::Type{Distribution{V, GaussianMeanPrecision}}, dist::Distribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType
    w = deepcopy(dist.params[:w])
    m = cholinv(w)*dist.params[:xi]

    return Distribution(V, GaussianMeanPrecision, m=m, w=w)
end

function convert(::Type{Distribution{V, GaussianMeanVariance}}, dist::Distribution{V, GaussianMeanPrecision}) where V<:VariateType
    v = cholinv(dist.params[:w])
    m = deepcopy(dist.params[:m])

    return Distribution(V, GaussianMeanVariance, m=m, v=v)
end

function convert(::Type{Distribution{V, GaussianMeanVariance}}, dist::Distribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType
    v = cholinv(dist.params[:w])
    m = v*dist.params[:xi]

    return Distribution(V, GaussianMeanVariance, m=m, v=v)
end

function convert(::Type{Distribution{V, GaussianWeightedMeanPrecision}}, dist::Distribution{V, GaussianMeanPrecision}) where V<:VariateType
    w = deepcopy(dist.params[:w])
    xi = w*dist.params[:m]

    return Distribution(V, GaussianWeightedMeanPrecision, xi=xi, w=w)
end

function convert(::Type{Distribution{V, GaussianWeightedMeanPrecision}}, dist::Distribution{V, GaussianMeanVariance}) where V<:VariateType
    w = cholinv(dist.params[:v])
    xi = w*dist.params[:m]

    return Distribution(V, GaussianWeightedMeanPrecision, xi=xi, w=w)
end

# Convert VariateTypes
convert(::Type{Distribution{Multivariate, GaussianMeanVariance}}, dist::Distribution{Univariate, GaussianMeanVariance}) =
    Distribution(Multivariate, GaussianMeanVariance, m=[dist.params[:m]], v=mat(dist.params[:v]))
convert(::Type{Distribution{Multivariate, GaussianMeanPrecision}}, dist::Distribution{Univariate, GaussianMeanPrecision}) =
    Distribution(Multivariate, GaussianMeanPrecision, m=[dist.params[:m]], w=mat(dist.params[:w]))
convert(::Type{Distribution{Multivariate, GaussianWeightedMeanPrecision}}, dist::Distribution{Univariate, GaussianWeightedMeanPrecision}) =
    Distribution(Multivariate, GaussianWeightedMeanPrecision, xi=[dist.params[:xi]], w=mat(dist.params[:w]))

function prod!(
    x::Distribution{Univariate, <:Gaussian},
    y::Distribution{Univariate, <:Gaussian},
    z::Distribution{Univariate, GaussianWeightedMeanPrecision}=Distribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=1.0))

    z.params[:xi] = unsafeWeightedMean(x) + unsafeWeightedMean(y)
    z.params[:w] = unsafePrecision(x) + unsafePrecision(y)

    return z
end

@symmetrical function prod!(
    x::Distribution{Univariate, <:Gaussian},
    y::Distribution{Univariate, PointMass},
    z::Distribution{Univariate, PointMass}=Distribution(Univariate, PointMass, m=0.0))

    z.params[:m] = y.params[:m]
    return z
end

function prod!(
    x::Distribution{Multivariate, <:Gaussian},
    y::Distribution{Multivariate, <:Gaussian},
    z::Distribution{Multivariate, GaussianWeightedMeanPrecision}=Distribution(Multivariate, GaussianWeightedMeanPrecision, xi=[NaN], w=transpose([NaN])))

    z.params[:xi] = unsafeWeightedMean(x) + unsafeWeightedMean(y)
    z.params[:w] = unsafePrecision(x) + unsafePrecision(y)

    return z
end

@symmetrical function prod!(
    x::Distribution{Multivariate, <:Gaussian},
    y::Distribution{Multivariate, PointMass},
    z::Distribution{Multivariate, PointMass}=Distribution(Multivariate, PointMass, m=[NaN]))

    z.params[:m] = deepcopy(y.params[:m])

    return z
end

function sample(dist::Distribution{Univariate, <:Gaussian})
    isProper(dist) || error("Cannot sample from improper distribution")
    (m,v) = unsafeMeanCov(dist)
    return sqrt(v)*randn() + m
end

function sample(dist::Distribution{Univariate, <:Gaussian}, n_samples::Int64)
    isProper(dist) || error("Cannot sample from improper distribution")
    (m,v) = unsafeMeanCov(dist)

    return sqrt(v).*randn(n_samples) .+ m
end

function sample(dist::Distribution{Multivariate, <:Gaussian})
    isProper(dist) || error("Cannot sample from improper distribution")
    (m,V) = unsafeMeanCov(dist)
    return (cholesky(default_cholesky_mode, V)).U' *randn(dims(dist)) + m
end

function sample(dist::Distribution{Multivariate, <:Gaussian}, n_samples::Int64)
    isProper(dist) || error("Cannot sample from improper distribution")
    (m,V) = unsafeMeanCov(dist)
    U = (cholesky(default_cholesky_mode, V)).U
    d = dims(dist)

    return [U' *randn(d) + m for i in 1:n_samples]
end

function naturalParams(dist::Distribution{<:VariateType, <:Gaussian})
    (xi, w) = unsafeWeightedMeanPrecision(dist)
    return vcat(xi, -0.5*vec(w))
end

standardDistribution(::Type{Univariate}, ::Type{<:Gaussian}; η::Vector) = Distribution(Univariate, GaussianWeightedMeanPrecision, xi=η[1], w=-2*η[2])
function standardDistribution(::Type{Multivariate}, ::Type{<:Gaussian}; η::Vector)
    d = Int(-0.5 + 0.5*sqrt(1 + 4*length(η))) # Extract dimensionality
    η_1 = η[1:d]
    η_2 = reshape(η[d+1:end], d, d)
    return Distribution(Multivariate, GaussianWeightedMeanPrecision, xi=η_1, w=-2.0*η_2)
end

logNormalizer(::Type{Univariate}, ::Type{<:Gaussian}; η::Vector) = -η[1]^2/(4*η[2]) - 0.5*log(-2*η[2])
function logNormalizer(::Type{Multivariate}, ::Type{<:Gaussian}; η::Vector) 
    d = Int(-0.5 + 0.5*sqrt(1 + 4*length(η))) # Extract dimensionality
    η_1 = η[1:d]
    η_2 = reshape(η[d+1:end], d, d)
    return η_1'*pinv(-4*η_2)*η_1 - 0.5*logdet(-2*η_2)
end

logPdf(V::Type{Univariate}, ::Type{F}, x::Number; η::Vector) where F<:Gaussian = -0.5*log(2pi) + vcat(x, x^2)'*η - logNormalizer(V, F; η=η)
function logPdf(V::Type{Multivariate}, ::Type{F}, x::Vector; η::Vector) where F<:Gaussian
    d = Int(-0.5 + 0.5*sqrt(1 + 4*length(η))) # Extract dimensionality
    return -0.5*d*log(2pi) + vcat(x, vec(x*x'))'*η - logNormalizer(V, F; η=η)
end

# Entropy functional
function differentialEntropy(dist::Distribution{Univariate, <:Gaussian})
    return  0.5*log(unsafeCov(dist)) +
            0.5*log(2*pi) +
            0.5
end

function differentialEntropy(dist::Distribution{Multivariate, <:Gaussian})
    d = dims(dist)[1]
    return  0.5*logdet(unsafeCov(dist)) +
            (d/2)*log(2*pi) +
            d/2
end
