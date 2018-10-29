export GaussianWeightedMeanPrecision

abstract type GaussianWeightedMeanPrecision <: Gaussian end

slug(::Type{GaussianWeightedMeanPrecision}) = "𝒩"

format(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = "$(slug(GaussianWeightedMeanPrecision))(xi=$(format(dist.params[:xi])), w=$(format(dist.params[:w])))"

ProbabilityDistribution(::Type{Univariate}, ::Type{GaussianWeightedMeanPrecision}; xi=0.0, w=1.0) = ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}(Dict(:xi=>xi, :w=>w))
ProbabilityDistribution(::Type{GaussianWeightedMeanPrecision}; xi::Number=0.0, w::Number=1.0) = ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}(Dict(:xi=>xi, :w=>w))
ProbabilityDistribution(::Type{Multivariate}, ::Type{GaussianWeightedMeanPrecision}; xi=[0.0], w=transpose([1.0])) = ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}(Dict(:xi=>xi, :w=>w))

dims(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = length(dist.params[:xi])

vague(::Type{GaussianWeightedMeanPrecision}) = ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=tiny)
vague(::Type{GaussianWeightedMeanPrecision}, dims::Int64) = ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=zeros(dims), w=tiny*diageye(dims))

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