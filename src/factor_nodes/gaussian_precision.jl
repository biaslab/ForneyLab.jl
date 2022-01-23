format(dist::Distribution{V, GaussianMeanPrecision}) where V<:VariateType = "$(slug(GaussianMeanPrecision))(m=$(format(dist.params[:m])), w=$(format(dist.params[:w])))"

Distribution(::Type{Univariate}, ::Type{GaussianMeanPrecision}; m=0.0, w=1.0) = Distribution{Univariate, GaussianMeanPrecision}(Dict(:m=>m, :w=>w))
Distribution(::Type{GaussianMeanPrecision}; m::Number=0.0, w::Number=1.0) = Distribution{Univariate, GaussianMeanPrecision}(Dict(:m=>m, :w=>w))
Distribution(::Type{Multivariate}, ::Type{GaussianMeanPrecision}; m=[0.0], w=transpose([1.0])) = Distribution{Multivariate, GaussianMeanPrecision}(Dict(:m=>m, :w=>w))

dims(dist::Distribution{V, GaussianMeanPrecision}) where V<:VariateType = size(dist.params[:m])

vague(::Type{GaussianMeanPrecision}) = Distribution(Univariate, GaussianMeanPrecision, m=0.0, w=tiny)
vague(::Type{GaussianMeanPrecision}, dims::Tuple{Int64}) = Distribution(Multivariate, GaussianMeanPrecision, m=zeros(dims), w=tiny*diageye(dims[1]))

unsafeMean(dist::Distribution{V, GaussianMeanPrecision}) where V<:VariateType = deepcopy(dist.params[:m]) # unsafe mean

unsafeMode(dist::Distribution{V, GaussianMeanPrecision}) where V<:VariateType = deepcopy(dist.params[:m]) # unsafe mode

unsafeVar(dist::Distribution{Univariate, GaussianMeanPrecision}) = 1.0/dist.params[:w] # unsafe variance
unsafeVar(dist::Distribution{Multivariate, GaussianMeanPrecision}) = diag(cholinv(dist.params[:w]))

unsafeCov(dist::Distribution{V, GaussianMeanPrecision}) where V<:VariateType = cholinv(dist.params[:w])

unsafeMeanCov(dist::Distribution{V, GaussianMeanPrecision}) where V<:VariateType = (deepcopy(dist.params[:m]), cholinv(dist.params[:w]))

unsafeWeightedMean(dist::Distribution{V, GaussianMeanPrecision}) where V<:VariateType = dist.params[:w]*dist.params[:m] # unsafe weighted mean

unsafePrecision(dist::Distribution{V, GaussianMeanPrecision}) where V<:VariateType = deepcopy(dist.params[:w]) # unsafe precision

unsafeMeanPrecision(dist::Distribution{V, GaussianMeanPrecision}) where V<:VariateType = (deepcopy(dist.params[:m]), deepcopy(dist.params[:w]))

unsafeWeightedMeanPrecision(dist::Distribution{V, GaussianMeanPrecision}) where V<:VariateType = (dist.params[:w]*dist.params[:m], deepcopy(dist.params[:w]))

logPdf(dist::Distribution{Univariate, GaussianMeanPrecision}, x) = -0.5*(log(2pi) - log(dist.params[:w]) + (x-dist.params[:m])^2*dist.params[:w])
logPdf(dist::Distribution{Multivariate, GaussianMeanPrecision}, x) = -0.5*(dims(dist)[1]*log(2pi) - logdet(dist.params[:w]) + transpose(x-dist.params[:m])*dist.params[:w]*(x-dist.params[:m]))

isProper(dist::Distribution{Univariate, GaussianMeanPrecision}) = (floatmin(Float64) < dist.params[:w] < floatmax(Float64))
isProper(dist::Distribution{Multivariate, GaussianMeanPrecision}) = isRoundedPosDef(dist.params[:w])

function ==(t::Distribution{V, GaussianMeanPrecision}, u::Distribution{V, GaussianMeanPrecision}) where V<:VariateType
    (t === u) && return true
    isApproxEqual(t.params[:m], u.params[:m]) && isApproxEqual(t.params[:w], u.params[:w])
end

function ==(t::Distribution{V, GaussianMeanPrecision}, u::Distribution{V, F}) where {V<:VariateType, F<:Gaussian}
    (t === u) && return true
    isApproxEqual(t.params[:m], unsafeMean(u)) && isApproxEqual(t.params[:w], unsafePrecision(u))
end

# Average energy functional
function averageEnergy(::Type{GaussianMeanPrecision}, marg_out::Distribution{Univariate}, marg_mean::Distribution{Univariate}, marg_prec::Distribution{Univariate})
    (m_mean, v_mean) = unsafeMeanCov(marg_mean)
    (m_out, v_out) = unsafeMeanCov(marg_out)

    0.5*log(2*pi) -
    0.5*unsafeLogMean(marg_prec) +
    0.5*unsafeMean(marg_prec)*(v_out + v_mean + (m_out - m_mean)^2)
end

function averageEnergy(::Type{GaussianMeanPrecision}, marg_out::Distribution{Multivariate}, marg_mean::Distribution{Multivariate}, marg_prec::Distribution{MatrixVariate})
    (m_mean, v_mean) = unsafeMeanCov(marg_mean)
    (m_out, v_out) = unsafeMeanCov(marg_out)

    0.5*dims(marg_out)[1]*log(2*pi) -
    0.5*unsafeDetLogMean(marg_prec) +
    0.5*tr( unsafeMean(marg_prec)*(v_out + v_mean + (m_out - m_mean)*(m_out - m_mean)' ))
end

function averageEnergy(::Type{GaussianMeanPrecision}, marg_out_mean::Distribution{Multivariate, <:Gaussian}, marg_prec::Distribution{Univariate})
    (m, V) = unsafeMeanCov(marg_out_mean)

    0.5*log(2*pi) -
    0.5*unsafeLogMean(marg_prec) +
    0.5*unsafeMean(marg_prec)*( V[1,1] - V[1,2] - V[2,1] + V[2,2] + (m[1] - m[2])^2 )
end

function averageEnergy(::Type{GaussianMeanPrecision}, marg_out_mean::Distribution{Multivariate, <:Gaussian}, marg_prec::Distribution{MatrixVariate})
    (m, V) = unsafeMeanCov(marg_out_mean)
    d = Int64(dims(marg_out_mean)[1]/2)

    0.5*d*log(2*pi) -
    0.5*unsafeDetLogMean(marg_prec) +
    0.5*tr( unsafeMean(marg_prec)*( V[1:d,1:d] - V[1:d,d+1:end] - V[d+1:end,1:d] + V[d+1:end,d+1:end] + (m[1:d] - m[d+1:end])*(m[1:d] - m[d+1:end])' ) )
end
