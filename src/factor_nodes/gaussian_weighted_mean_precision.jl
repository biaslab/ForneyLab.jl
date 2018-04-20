export GaussianWeightedMeanPrecision

abstract type GaussianWeightedMeanPrecision <: Gaussian end

slug(::Type{GaussianWeightedMeanPrecision}) = "𝒩"

format{V<:VariateType}(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) = "$(slug(GaussianWeightedMeanPrecision))(m=$(format(dist.params[:m])), xi=$(format(dist.params[:xi])))"

ProbabilityDistribution(::Type{Univariate}, ::Type{GaussianWeightedMeanPrecision}; xi=0.0, w=1.0) = ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}(Dict(xi=xi, w=w))
ProbabilityDistribution(::Type{GaussianWeightedMeanPrecision}; xi::Number=0.0, w::Number=1.0) = ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}(Dict(xi=xi, w=w))
ProbabilityDistribution(::Type{Multivariate}, ::Type{GaussianWeightedMeanPrecision}; xi=[0.0], w=[1.0].') = ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}(Dict(xi=xi, w=w))

dims{V<:VariateType}(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) = length(dist.params[:m])

vague(::Type{GaussianWeightedMeanPrecision}) = ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=tiny)
vague(::Type{GaussianWeightedMeanPrecision}, dims::Int64) = ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=zeros(dims), w=tiny*diageye(dims))

unsafeMean(dist::ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}) = dist.params[:xi]/dist.params[:w] # unsafe mean
unsafeMean(dist::ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}) = cholinv(dist.params[:w])*dist.params[:xi]

unsafeVar(dist::ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}) = 1.0/dist.params[:w] # unsafe variance
unsafeVar(dist::ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}) = diag(cholinv(dist.params[:w]))

unsafeCov(dist::ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}) = 1.0/dist.params[:w] # unsafe covariance
unsafeCov(dist::ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}) = cholinv(dist.params[:w])

unsafeWeightedMean{V<:VariateType}(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) = deepcopy(dist.params[:xi]) # unsafe weighted mean

unsafePrecision{V<:VariateType}(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) = deepcopy(dist.params[:w]) # unsafe precision

isProper(dist::ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}) = (realmin(Float64) < dist.params[:w] < realmax(Float64))
isProper(dist::ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}) = isRoundedPosDef(dist.params[:w])

function =={V<:VariateType}(t::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}, u::ProbabilityDistribution{V, GaussianWeightedMeanPrecision})
    (t === u) && return true
    isApproxEqual(t.params[:xi], u.params[:xi]) && isApproxEqual(t.params[:w], u.params[:w])
end

function =={V<:VariateTypem F<:Gaussian}(t::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}, u::ProbabilityDistribution{V, F})
    (t === u) && return true
    isApproxEqual(t.params[:xi], unsafeWeightedMean(u)) && isApproxEqual(t.params[:w], unsafePrecision(u))
end