format(dist::Distribution{V, Gaussian{Canonical}}) where V<:VariateType = "$(slug(Gaussian{Canonical}))(xi=$(format(dist.params[:xi])), w=$(format(dist.params[:w])))"

Distribution(::Type{Univariate}, ::Type{Gaussian{Canonical}}; xi=0.0, w=1.0) = Distribution{Univariate, Gaussian{Canonical}}(Dict(:xi=>xi, :w=>w))
Distribution(::Type{Gaussian{Canonical}}; xi::Number=0.0, w::Number=1.0) = Distribution{Univariate, Gaussian{Canonical}}(Dict(:xi=>xi, :w=>w))
Distribution(::Type{Multivariate}, ::Type{Gaussian{Canonical}}; xi=[0.0], w=transpose([1.0])) = Distribution{Multivariate, Gaussian{Canonical}}(Dict(:xi=>xi, :w=>w))

dims(dist::Distribution{V, Gaussian{Canonical}}) where V<:VariateType = size(dist.params[:xi])

vague(::Type{Gaussian{Canonical}}) = Distribution(Univariate, Gaussian{Canonical}, xi=0.0, w=tiny)
vague(::Type{Gaussian{Canonical}}, dims::Tuple{Int64}) = Distribution(Multivariate, Gaussian{Canonical}, xi=zeros(dims), w=tiny*diageye(dims[1]))

unsafeMean(dist::Distribution{V, Gaussian{Canonical}}) where V<:VariateType = cholinv(dist.params[:w])*dist.params[:xi]

unsafeMode(dist::Distribution{V, Gaussian{Canonical}}) where V<:VariateType = cholinv(dist.params[:w])*dist.params[:xi]

unsafeVar(dist::Distribution{Univariate, Gaussian{Canonical}}) = 1.0/dist.params[:w] # unsafe variance
unsafeVar(dist::Distribution{Multivariate, Gaussian{Canonical}}) = diag(cholinv(dist.params[:w]))

unsafeCov(dist::Distribution{V, Gaussian{Canonical}}) where V<:VariateType = cholinv(dist.params[:w])

function unsafeMeanCov(dist::Distribution{V, Gaussian{Canonical}}) where V<:VariateType
    v = cholinv(dist.params[:w])
    return (v*dist.params[:xi], v)
end

unsafeWeightedMean(dist::Distribution{V, Gaussian{Canonical}}) where V<:VariateType = deepcopy(dist.params[:xi]) # unsafe weighted mean

unsafePrecision(dist::Distribution{V, Gaussian{Canonical}}) where V<:VariateType = deepcopy(dist.params[:w]) # unsafe precision

function unsafeMeanPrecision(dist::Distribution{V, Gaussian{Canonical}}) where V<:VariateType
    v = cholinv(dist.params[:w])
    return (v*dist.params[:xi], dist.params[:w])
end

unsafeWeightedMeanPrecision(dist::Distribution{V, Gaussian{Canonical}}) where V<:VariateType = (deepcopy(dist.params[:xi]), deepcopy(dist.params[:w]))

function logPdf(dist::Distribution{Univariate, Gaussian{Canonical}}, x)
    m = dist.params[:xi]/dist.params[:w]
    return -0.5*(log(2pi) - log(dist.params[:w]) + (x-m)^2*dist.params[:w])
end

function logPdf(dist::Distribution{Multivariate, Gaussian{Canonical}}, x) 
    m = cholinv(dist.params[:w])*dist.params[:xi]
    return -0.5*(dims(dist)[1]*log(2pi) - logdet(dist.params[:w]) + transpose(x-m)*dist.params[:w]*(x-m))
end

isProper(dist::Distribution{Univariate, Gaussian{Canonical}}) = (floatmin(Float64) < dist.params[:w] < floatmax(Float64))
isProper(dist::Distribution{Multivariate, Gaussian{Canonical}}) = isRoundedPosDef(dist.params[:w])

function ==(t::Distribution{V, Gaussian{Canonical}}, u::Distribution{V, Gaussian{Canonical}}) where V<:VariateType
    (t === u) && return true
    isApproxEqual(t.params[:xi], u.params[:xi]) && isApproxEqual(t.params[:w], u.params[:w])
end

function ==(t::Distribution{V, Gaussian{Canonical}}, u::Distribution{V, <:Gaussian}) where V<:VariateType
    (t === u) && return true
    isApproxEqual(t.params[:xi], unsafeWeightedMean(u)) && isApproxEqual(t.params[:w], unsafePrecision(u))
end