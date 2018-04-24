export GaussianWeightedMeanPrecision

abstract type GaussianWeightedMeanPrecision <: Gaussian end

slug(::Type{GaussianWeightedMeanPrecision}) = "𝒩"

format{V<:VariateType}(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) = "$(slug(GaussianWeightedMeanPrecision))(xi=$(format(dist.params[:xi])), w=$(format(dist.params[:w])))"

ProbabilityDistribution(::Type{Univariate}, ::Type{GaussianWeightedMeanPrecision}; xi=0.0, w=1.0) = ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}(Dict(:xi=>xi, :w=>w))
ProbabilityDistribution(::Type{GaussianWeightedMeanPrecision}; xi::Number=0.0, w::Number=1.0) = ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}(Dict(:xi=>xi, :w=>w))
ProbabilityDistribution(::Type{Multivariate}, ::Type{GaussianWeightedMeanPrecision}; xi=[0.0], w=[1.0].') = ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}(Dict(:xi=>xi, :w=>w))

dims{V<:VariateType}(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) = length(dist.params[:xi])

vague(::Type{GaussianWeightedMeanPrecision}) = ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=tiny)
vague(::Type{GaussianWeightedMeanPrecision}, dims::Int64) = ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=zeros(dims), w=tiny*diageye(dims))

unsafeMean{V<:VariateType}(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) = cholinv(dist.params[:w])*dist.params[:xi]

unsafeVar(dist::ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}) = 1.0/dist.params[:w] # unsafe variance
unsafeVar(dist::ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}) = diag(cholinv(dist.params[:w]))

unsafeCov{V<:VariateType}(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) = cholinv(dist.params[:w])

function unsafeMeanCov{V<:VariateType}(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision})
    v = cholinv(dist.params[:w])
    return (v*dist.params[:xi], v)
end

unsafeWeightedMean{V<:VariateType}(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) = deepcopy(dist.params[:xi]) # unsafe weighted mean

unsafePrecision{V<:VariateType}(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) = deepcopy(dist.params[:w]) # unsafe precision

unsafeWeightedMeanPrecision{V<:VariateType}(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) = (deepcopy(dist.params[:xi]), deepcopy(dist.params[:w]))

isProper(dist::ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}) = (realmin(Float64) < dist.params[:w] < realmax(Float64))
isProper(dist::ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}) = isRoundedPosDef(dist.params[:w])

function =={V<:VariateType}(t::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}, u::ProbabilityDistribution{V, GaussianWeightedMeanPrecision})
    (t === u) && return true
    isApproxEqual(t.params[:xi], u.params[:xi]) && isApproxEqual(t.params[:w], u.params[:w])
end

function =={V<:VariateType, F<:Gaussian}(t::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}, u::ProbabilityDistribution{V, F})
    (t === u) && return true
    isApproxEqual(t.params[:xi], unsafeWeightedMean(u)) && isApproxEqual(t.params[:w], unsafePrecision(u))
end