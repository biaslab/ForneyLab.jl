export GaussianWeightedMeanPrecision

"""
Description:

    A Gaussian with weighted-mean-precision parameterization:

    f(out,xi,w) = 𝒩(out|xi,w) = (2π)^{-D/2} |w|^{1/2} exp(-1/2 xi' w^{-1} xi) exp(-1/2 xi' w xi + out' xi)

Interfaces:

    1. out
    2. xi (weighted mean, w*m)
    3. w (precision)

Construction:

    GaussianWeightedMeanPrecision(out, xi, w, id=:some_id)
"""
mutable struct GaussianWeightedMeanPrecision <: Gaussian
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function GaussianWeightedMeanPrecision(out, xi, w; id=generateId(GaussianWeightedMeanPrecision))
        @ensureVariables(out, xi, w)
        self = new(id, Array{Interface}(undef, 3), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:xi] = self.interfaces[2] = associate!(Interface(self), xi)
        self.i[:w] = self.interfaces[3] = associate!(Interface(self), w)

        return self
    end
end

slug(::Type{GaussianWeightedMeanPrecision}) = "𝒩"

format(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = "$(slug(GaussianWeightedMeanPrecision))(xi=$(format(dist.params[:xi])), w=$(format(dist.params[:w])))"

ProbabilityDistribution(::Type{Univariate}, ::Type{GaussianWeightedMeanPrecision}; xi=0.0, w=1.0) = ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}(Dict(:xi=>xi, :w=>w))
ProbabilityDistribution(::Type{GaussianWeightedMeanPrecision}; xi::Number=0.0, w::Number=1.0) = ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}(Dict(:xi=>xi, :w=>w))
ProbabilityDistribution(::Type{Multivariate}, ::Type{GaussianWeightedMeanPrecision}; xi=[0.0], w=transpose([1.0])) = ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}(Dict(:xi=>xi, :w=>w))

dims(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = length(dist.params[:xi])

vague(::Type{GaussianWeightedMeanPrecision}) = ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=tiny)
vague(::Type{GaussianWeightedMeanPrecision}, dims::Int64) = ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=zeros(dims), w=tiny*diageye(dims))
vague(::Type{GaussianWeightedMeanPrecision}, dims::Tuple{Int64}) = ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=zeros(dims), w=tiny*diageye(dims[1]))

unsafeMean(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = cholinv(dist.params[:w])*dist.params[:xi]

unsafeMode(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = cholinv(dist.params[:w])*dist.params[:xi]

unsafeVar(dist::ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}) = 1.0/dist.params[:w] # unsafe variance
unsafeVar(dist::ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}) = diag(cholinv(dist.params[:w]))

unsafeCov(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = cholinv(dist.params[:w])

function unsafeMeanCov(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType
    v = cholinv(dist.params[:w])
    return (v*dist.params[:xi], v)
end

unsafeWeightedMean(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = deepcopy(dist.params[:xi]) # unsafe weighted mean

unsafePrecision(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = deepcopy(dist.params[:w]) # unsafe precision

unsafeWeightedMeanPrecision(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = (deepcopy(dist.params[:xi]), deepcopy(dist.params[:w]))

function logPdf(dist::ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}, x)
    m = dist.params[:xi]/dist.params[:w]
    return -0.5*(log(2pi) - log(dist.params[:w]) + (x-m)^2*dist.params[:w])
end

function logPdf(dist::ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}, x) 
    m = cholinv(dist.params[:w])*dist.params[:xi]
    return -0.5*(dims(dist)*log(2pi) - log(det(dist.params[:w])) + transpose(x-m)*dist.params[:w]*(x-m))
end

isProper(dist::ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}) = (floatmin(Float64) < dist.params[:w] < floatmax(Float64))
isProper(dist::ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}) = isRoundedPosDef(dist.params[:w])

function ==(t::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}, u::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType
    (t === u) && return true
    isApproxEqual(t.params[:xi], u.params[:xi]) && isApproxEqual(t.params[:w], u.params[:w])
end

function ==(t::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}, u::ProbabilityDistribution{V, F}) where {V<:VariateType, F<:Gaussian}
    (t === u) && return true
    isApproxEqual(t.params[:xi], unsafeWeightedMean(u)) && isApproxEqual(t.params[:w], unsafePrecision(u))
end