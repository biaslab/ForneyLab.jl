format(dist::Distribution{V, Gaussian{Precision}}) where V<:VariateType = "$(slug(Gaussian{Precision}))(m=$(format(dist.params[:m])), w=$(format(dist.params[:w])))"

Distribution(::Type{Univariate}, ::Type{Gaussian{Precision}}; m=0.0, w=1.0) = Distribution{Univariate, Gaussian{Precision}}(Dict(:m=>m, :w=>w))
Distribution(::Type{Gaussian{Precision}}; m::Number=0.0, w::Number=1.0) = Distribution{Univariate, Gaussian{Precision}}(Dict(:m=>m, :w=>w))
Distribution(::Type{Multivariate}, ::Type{Gaussian{Precision}}; m=[0.0], w=transpose([1.0])) = Distribution{Multivariate, Gaussian{Precision}}(Dict(:m=>m, :w=>w))

dims(dist::Distribution{V, Gaussian{Precision}}) where V<:VariateType = size(dist.params[:m])

vague(::Type{Gaussian{Precision}}) = Distribution(Univariate, Gaussian{Precision}, m=0.0, w=tiny)
vague(::Type{Gaussian{Precision}}, dims::Tuple{Int64}) = Distribution(Multivariate, Gaussian{Precision}, m=zeros(dims), w=tiny*diageye(dims[1]))

unsafeMean(dist::Distribution{V, Gaussian{Precision}}) where V<:VariateType = deepcopy(dist.params[:m]) # unsafe mean

unsafeMode(dist::Distribution{V, Gaussian{Precision}}) where V<:VariateType = deepcopy(dist.params[:m]) # unsafe mode

unsafeVar(dist::Distribution{Univariate, Gaussian{Precision}}) = 1.0/dist.params[:w] # unsafe variance
unsafeVar(dist::Distribution{Multivariate, Gaussian{Precision}}) = diag(cholinv(dist.params[:w]))

unsafeCov(dist::Distribution{V, Gaussian{Precision}}) where V<:VariateType = cholinv(dist.params[:w])

unsafeMeanCov(dist::Distribution{V, Gaussian{Precision}}) where V<:VariateType = (deepcopy(dist.params[:m]), cholinv(dist.params[:w]))

unsafeWeightedMean(dist::Distribution{V, Gaussian{Precision}}) where V<:VariateType = dist.params[:w]*dist.params[:m] # unsafe weighted mean

unsafePrecision(dist::Distribution{V, Gaussian{Precision}}) where V<:VariateType = deepcopy(dist.params[:w]) # unsafe precision

unsafeMeanPrecision(dist::Distribution{V, Gaussian{Precision}}) where V<:VariateType = (deepcopy(dist.params[:m]), deepcopy(dist.params[:w]))

unsafeWeightedMeanPrecision(dist::Distribution{V, Gaussian{Precision}}) where V<:VariateType = (dist.params[:w]*dist.params[:m], deepcopy(dist.params[:w]))

logPdf(dist::Distribution{Univariate, Gaussian{Precision}}, x) = -0.5*(log(2pi) - log(dist.params[:w]) + (x-dist.params[:m])^2*dist.params[:w])
logPdf(dist::Distribution{Multivariate, Gaussian{Precision}}, x) = -0.5*(dims(dist)[1]*log(2pi) - logdet(dist.params[:w]) + transpose(x-dist.params[:m])*dist.params[:w]*(x-dist.params[:m]))

isProper(dist::Distribution{Univariate, Gaussian{Precision}}) = (floatmin(Float64) < dist.params[:w] < floatmax(Float64))
isProper(dist::Distribution{Multivariate, Gaussian{Precision}}) = isRoundedPosDef(dist.params[:w])

function ==(t::Distribution{V, Gaussian{Precision}}, u::Distribution{V, Gaussian{Precision}}) where V<:VariateType
    (t === u) && return true
    isApproxEqual(t.params[:m], u.params[:m]) && isApproxEqual(t.params[:w], u.params[:w])
end

function ==(t::Distribution{V, Gaussian{Precision}}, u::Distribution{V, F}) where {V<:VariateType, F<:Gaussian}
    (t === u) && return true
    isApproxEqual(t.params[:m], unsafeMean(u)) && isApproxEqual(t.params[:w], unsafePrecision(u))
end

# Average energy functional
function averageEnergy(::Type{Gaussian{Precision}}, marg_out::Distribution{Univariate}, marg_mean::Distribution{Univariate}, marg_prec::Distribution{Univariate})
    (m_mean, v_mean) = unsafeMeanCov(marg_mean)
    (m_out, v_out) = unsafeMeanCov(marg_out)

    0.5*log(2*pi) -
    0.5*unsafeLogMean(marg_prec) +
    0.5*unsafeMean(marg_prec)*(v_out + v_mean + (m_out - m_mean)^2)
end

function averageEnergy(::Type{Gaussian{Precision}}, marg_out::Distribution{Multivariate}, marg_mean::Distribution{Multivariate}, marg_prec::Distribution{MatrixVariate})
    (m_mean, v_mean) = unsafeMeanCov(marg_mean)
    (m_out, v_out) = unsafeMeanCov(marg_out)

    0.5*dims(marg_out)[1]*log(2*pi) -
    0.5*unsafeDetLogMean(marg_prec) +
    0.5*tr( unsafeMean(marg_prec)*(v_out + v_mean + (m_out - m_mean)*(m_out - m_mean)' ))
end

function averageEnergy(::Type{Gaussian{Precision}}, marg_out_mean::Distribution{Multivariate, <:Gaussian}, marg_prec::Distribution{Univariate})
    (m, V) = unsafeMeanCov(marg_out_mean)

    0.5*log(2*pi) -
    0.5*unsafeLogMean(marg_prec) +
    0.5*unsafeMean(marg_prec)*( V[1,1] - V[1,2] - V[2,1] + V[2,2] + (m[1] - m[2])^2 )
end

function averageEnergy(::Type{Gaussian{Precision}}, marg_out_mean::Distribution{Multivariate, <:Gaussian}, marg_prec::Distribution{MatrixVariate})
    (m, V) = unsafeMeanCov(marg_out_mean)
    d = Int64(dims(marg_out_mean)[1]/2)

    0.5*d*log(2*pi) -
    0.5*unsafeDetLogMean(marg_prec) +
    0.5*tr( unsafeMean(marg_prec)*( V[1:d,1:d] - V[1:d,d+1:end] - V[d+1:end,1:d] + V[d+1:end,d+1:end] + (m[1:d] - m[d+1:end])*(m[1:d] - m[d+1:end])' ) )
end
