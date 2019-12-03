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

function format(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType
    return "$(slug(GaussianWeightedMeanPrecision))(xi=$(format(dist.params[:xi])), w=$(format(dist.params[:w])))"
end

function ProbabilityDistribution(::Type{Univariate}, ::Type{GaussianWeightedMeanPrecision}; xi=0.0, w=1.0)
    return ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}(Dict(:xi=>xi, :w=>w))
end

function ProbabilityDistribution(::Type{GaussianWeightedMeanPrecision}; xi::Number=0.0, w::Number=1.0)
    return ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}(Dict(:xi=>xi, :w=>w))
end

function ProbabilityDistribution(::Type{Multivariate}, ::Type{GaussianWeightedMeanPrecision}; xi=[0.0], w=Matrix{Float64}(I,1,1))
    if isa(w, PDMat)
        return ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}(Dict(:xi=>xi, :w=>w))
    else
        return ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}(Dict(:xi=>xi, :w=>PDMat(Matrix(w))))
    end
end

function dims(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType
    return length(dist.params[:xi])
end

function vague(::Type{GaussianWeightedMeanPrecision})
    return ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=tiny)
end

function vague(::Type{GaussianWeightedMeanPrecision}, dims::Int64)
    return ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=zeros(dims), w=PDMat(tiny*Matrix{Float64}(I,dims,dims)))
end

function vague(::Type{GaussianWeightedMeanPrecision}, dims::Tuple{Int64})
    return ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=zeros(dims), w=PDMat(tiny*Matrix{Float64}(I,dims[1],dims[1])))
end

function unsafeMode(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType
    return inv(dist.params[:w])*dist.params[:xi]
end

function unsafeMean(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType
    return inv(dist.params[:w])*dist.params[:xi]
end

function unsafeVar(dist::ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision})
    return 1.0/dist.params[:w]
end

function unsafeVar(dist::ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision})
    return diag(inv(dist.params[:w]))
end

function unsafePrecision(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType
    return deepcopy(dist.params[:w])
end

function unsafeCov(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType
    return inv(dist.params[:w])
end

function unsafeMeanCov(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType
    v = inv(dist.params[:w])
    return (v*dist.params[:xi], v)
end

function unsafeWeightedMean(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType
    return deepcopy(dist.params[:xi])
end

function unsafeWeightedMeanPrecision(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType
    return (deepcopy(dist.params[:xi]), deepcopy(dist.params[:w]))
end

function isProper(dist::ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision})
    return (floatmin(Float64) < dist.params[:w] < floatmax(Float64))
end

function isProper(dist::ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision})
    return isRoundedPosDef(dist.params[:w])
end

function ==(t::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}, u::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType
    (t === u) && return true
    isApproxEqual(t.params[:xi], u.params[:xi]) && isApproxEqual(t.params[:w], u.params[:w])
end

function ==(t::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}, u::ProbabilityDistribution{V, F}) where {V<:VariateType, F<:Gaussian}
    (t === u) && return true
    isApproxEqual(t.params[:xi], unsafeWeightedMean(u)) && isApproxEqual(t.params[:w], unsafePrecision(u))
end
