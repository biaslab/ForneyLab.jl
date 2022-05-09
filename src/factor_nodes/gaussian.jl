export 
Gaussian, 
Moments, 
Precision, 
Canonical, 
prod!, 
convert,
GaussianMeanVariance,
GaussianMeanPrecision,
GaussianWeightedMeanPrecision

abstract type GaussianParameterization end
abstract type Moments <: GaussianParameterization end
abstract type Precision <: GaussianParameterization end
abstract type Canonical <: GaussianParameterization end

interfaceHandles(::Type{Moments}) = (:m, :v)
interfaceHandles(::Type{Precision}) = (:m, :w)
interfaceHandles(::Type{Canonical}) = (:xi, :w)

"""
Description:

    A Gaussian with moments, precision or canonical parameterization:

    f(out, m, v) = 𝒩(out | m, v)
    f(out, m, w) = 𝒩(out | m, w^{-1})
    f(out, xi, w) = 𝒩(out | w^{-1}*xi, w^{-1})

Interfaces:

    1. out
    2. m, xi
    3. v, w

Construction:

    Gaussian(out, m, v, id=:some_id)
    Gaussian{Moments}(out, m, v, id=:some_id)
    Gaussian{Precision}(out, m, w, id=:some_id)
    Gaussian{Canonical}(out, xi, w, id=:some_id)
"""
mutable struct Gaussian{T<:GaussianParameterization} <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    # Default to Moments parameterization
    function Gaussian(out, m, v; id=generateId(Gaussian{Moments}))
        @ensureVariables(out, m, v)
        self = new{Moments}(id, Array{Interface}(undef, 3), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:m] = self.interfaces[2] = associate!(Interface(self), m)
        self.i[:v] = self.interfaces[3] = associate!(Interface(self), v)

        return self
    end

    # User-defined parameterization
    function Gaussian{T}(out, in1, in2; id=generateId(Gaussian{T})) where T<:GaussianParameterization
        @ensureVariables(out, in1, in2)
        self = new(id, Array{Interface}(undef, 3), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        h = interfaceHandles(T) # Extract interface handle symbols from approximation type
        self.i[h[1]] = self.interfaces[2] = associate!(Interface(self), in1)
        self.i[h[2]] = self.interfaces[3] = associate!(Interface(self), in2)

        return self
    end
end

"""Aliases for Gaussian definitions"""
const GaussianMeanVariance = Gaussian{Moments}  # For backwards compatibility
const GaussianMeanPrecision = Gaussian{Precision}
const GaussianWeightedMeanPrecision = Gaussian{Canonical}

slug(::Type{<:Gaussian}) = "𝒩"

# Convert parameterizations
function convert(::Type{Distribution{V, Gaussian{Precision}}}, dist::Distribution{V, Gaussian{Moments}}) where V<:VariateType
    w = cholinv(dist.params[:v])
    m = deepcopy(dist.params[:m])

    return Distribution(V, Gaussian{Precision}, m=m, w=w)
end

function convert(::Type{Distribution{V, Gaussian{Precision}}}, dist::Distribution{V, Gaussian{Canonical}}) where V<:VariateType
    w = deepcopy(dist.params[:w])
    m = cholinv(w)*dist.params[:xi]

    return Distribution(V, Gaussian{Precision}, m=m, w=w)
end

function convert(::Type{Distribution{V, Gaussian{Moments}}}, dist::Distribution{V, Gaussian{Precision}}) where V<:VariateType
    v = cholinv(dist.params[:w])
    m = deepcopy(dist.params[:m])

    return Distribution(V, Gaussian{Moments}, m=m, v=v)
end

function convert(::Type{Distribution{V, Gaussian{Moments}}}, dist::Distribution{V, Gaussian{Canonical}}) where V<:VariateType
    v = cholinv(dist.params[:w])
    m = v*dist.params[:xi]

    return Distribution(V, Gaussian{Moments}, m=m, v=v)
end

function convert(::Type{Distribution{V, Gaussian{Canonical}}}, dist::Distribution{V, Gaussian{Precision}}) where V<:VariateType
    w = deepcopy(dist.params[:w])
    xi = w*dist.params[:m]

    return Distribution(V, Gaussian{Canonical}, xi=xi, w=w)
end

function convert(::Type{Distribution{V, Gaussian{Canonical}}}, dist::Distribution{V, Gaussian{Moments}}) where V<:VariateType
    w = cholinv(dist.params[:v])
    xi = w*dist.params[:m]

    return Distribution(V, Gaussian{Canonical}, xi=xi, w=w)
end

# Convert VariateTypes
convert(::Type{Distribution{Multivariate, Gaussian{Moments}}}, dist::Distribution{Univariate, Gaussian{Moments}}) =
    Distribution(Multivariate, Gaussian{Moments}, m=[dist.params[:m]], v=mat(dist.params[:v]))
convert(::Type{Distribution{Multivariate, Gaussian{Precision}}}, dist::Distribution{Univariate, Gaussian{Precision}}) =
    Distribution(Multivariate, Gaussian{Precision}, m=[dist.params[:m]], w=mat(dist.params[:w]))
convert(::Type{Distribution{Multivariate, Gaussian{Canonical}}}, dist::Distribution{Univariate, Gaussian{Canonical}}) =
    Distribution(Multivariate, Gaussian{Canonical}, xi=[dist.params[:xi]], w=mat(dist.params[:w]))

function prod!(
    x::Distribution{Univariate, <:Gaussian},
    y::Distribution{Univariate, <:Gaussian},
    z::Distribution{Univariate, Gaussian{Canonical}}=Distribution(Univariate, Gaussian{Canonical}, xi=0.0, w=1.0))

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
    z::Distribution{Multivariate, Gaussian{Canonical}}=Distribution(Multivariate, Gaussian{Canonical}, xi=[NaN], w=transpose([NaN])))

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

standardDistribution(::Type{Univariate}, ::Type{<:Gaussian}; η::Vector) = Distribution(Univariate, Gaussian{Canonical}, xi=η[1], w=-2*η[2])
function standardDistribution(::Type{Multivariate}, ::Type{<:Gaussian}; η::Vector)
    d = Int(-0.5 + 0.5*sqrt(1 + 4*length(η))) # Extract dimensionality
    η_1 = η[1:d]
    η_2 = reshape(η[d+1:end], d, d)
    return Distribution(Multivariate, Gaussian{Canonical}, xi=η_1, w=-2.0*η_2)
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
