export GaussianWeightedMeanPrecision

abstract type GaussianWeightedMeanPrecision <: Gaussian end

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

function ProbabilityDistribution(::Type{Multivariate}, ::Type{GaussianWeightedMeanPrecision}; xi=[0.0], w=diageye(1))
    return ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}(Dict(:xi=>xi, :w=>w))
end

function dims(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType
    return length(dist.params[:xi])
end

function vague(::Type{GaussianWeightedMeanPrecision})
    return ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=tiny)
end

function vague(::Type{GaussianWeightedMeanPrecision}, dims::Int64)
    return ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=zeros(dims), w=tiny*diageye(dims))
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

function unsafePrecision(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType
    return deepcopy(dist.params[:w])
end

function unsafeWeightedMeanPrecision(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType
    return (deepcopy(dist.params[:xi]), deepcopy(dist.params[:w]))
end

function isProper(dist::ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision})
    return (floatmin(Float64) < dist.params[:w] < floatmax(Float64))
end

# Method is obsolete, non-proper multivariate Gaussian distributions cannot be constructed
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
