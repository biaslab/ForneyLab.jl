format(dist::Distribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = "$(slug(GaussianWeightedMeanPrecision))(xi=$(format(dist.params[:xi])), w=$(format(dist.params[:w])))"

Distribution(::Type{Univariate}, ::Type{GaussianWeightedMeanPrecision}; xi=0.0, w=1.0) = Distribution{Univariate, GaussianWeightedMeanPrecision}(Dict(:xi=>xi, :w=>w))
Distribution(::Type{GaussianWeightedMeanPrecision}; xi::Number=0.0, w::Number=1.0) = Distribution{Univariate, GaussianWeightedMeanPrecision}(Dict(:xi=>xi, :w=>w))
Distribution(::Type{Multivariate}, ::Type{GaussianWeightedMeanPrecision}; xi=[0.0], w=transpose([1.0])) = Distribution{Multivariate, GaussianWeightedMeanPrecision}(Dict(:xi=>xi, :w=>w))

dims(dist::Distribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = size(dist.params[:xi])

vague(::Type{GaussianWeightedMeanPrecision}) = Distribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=tiny)
vague(::Type{GaussianWeightedMeanPrecision}, dims::Tuple{Int64}) = Distribution(Multivariate, GaussianWeightedMeanPrecision, xi=zeros(dims), w=tiny*diageye(dims[1]))

unsafeMean(dist::Distribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = cholinv(dist.params[:w])*dist.params[:xi]

unsafeMode(dist::Distribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = cholinv(dist.params[:w])*dist.params[:xi]

unsafeVar(dist::Distribution{Univariate, GaussianWeightedMeanPrecision}) = 1.0/dist.params[:w] # unsafe variance
unsafeVar(dist::Distribution{Multivariate, GaussianWeightedMeanPrecision}) = diag(cholinv(dist.params[:w]))

unsafeCov(dist::Distribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = cholinv(dist.params[:w])

function unsafeMeanCov(dist::Distribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType
    v = cholinv(dist.params[:w])
    return (v*dist.params[:xi], v)
end

unsafeWeightedMean(dist::Distribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = deepcopy(dist.params[:xi]) # unsafe weighted mean

unsafePrecision(dist::Distribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = deepcopy(dist.params[:w]) # unsafe precision

function unsafeMeanPrecision(dist::Distribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType
    v = cholinv(dist.params[:w])
    return (v*dist.params[:xi], dist.params[:w])
end

unsafeWeightedMeanPrecision(dist::Distribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = (deepcopy(dist.params[:xi]), deepcopy(dist.params[:w]))

function logPdf(dist::Distribution{Univariate, GaussianWeightedMeanPrecision}, x)
    m = dist.params[:xi]/dist.params[:w]
    return -0.5*(log(2pi) - log(dist.params[:w]) + (x-m)^2*dist.params[:w])
end

function logPdf(dist::Distribution{Multivariate, GaussianWeightedMeanPrecision}, x) 
    m = cholinv(dist.params[:w])*dist.params[:xi]
    return -0.5*(dims(dist)[1]*log(2pi) - logdet(dist.params[:w]) + transpose(x-m)*dist.params[:w]*(x-m))
end

isProper(dist::Distribution{Univariate, GaussianWeightedMeanPrecision}) = (floatmin(Float64) < dist.params[:w] < floatmax(Float64))
isProper(dist::Distribution{Multivariate, GaussianWeightedMeanPrecision}) = isRoundedPosDef(dist.params[:w])

function ==(t::Distribution{V, GaussianWeightedMeanPrecision}, u::Distribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType
    (t === u) && return true
    isApproxEqual(t.params[:xi], u.params[:xi]) && isApproxEqual(t.params[:w], u.params[:w])
end

function ==(t::Distribution{V, GaussianWeightedMeanPrecision}, u::Distribution{V, <:Gaussian}) where V<:VariateType
    (t === u) && return true
    isApproxEqual(t.params[:xi], unsafeWeightedMean(u)) && isApproxEqual(t.params[:w], unsafePrecision(u))
end