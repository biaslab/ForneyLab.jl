format(dist::Distribution{V, Gaussian{Moments}}) where V<:VariateType = "$(slug(Gaussian{Moments}))(m=$(format(dist.params[:m])), v=$(format(dist.params[:v])))"

Distribution(::Type{Univariate}, ::Type{Gaussian{Moments}}; m=0.0, v=1.0) = Distribution{Univariate, Gaussian{Moments}}(Dict(:m=>m, :v=>v))
Distribution(::Type{Gaussian{Moments}}; m::Number=0.0, v::Number=1.0) = Distribution{Univariate, Gaussian{Moments}}(Dict(:m=>m, :v=>v))
Distribution(::Type{Multivariate}, ::Type{Gaussian{Moments}}; m=[0.0], v=mat(1.0)) = Distribution{Multivariate, Gaussian{Moments}}(Dict(:m=>m, :v=>v))
# Default to Moments parameterization
Distribution(::Type{Univariate}, ::Type{Gaussian}; m=0.0, v=1.0) = Distribution(Univariate, Gaussian{Moments}, m=m, v=v)
Distribution(::Type{Gaussian}; m::Number=0.0, v::Number=1.0) = Distribution(Univariate, Gaussian{Moments}, m=m, v=v)
Distribution(::Type{Multivariate}, ::Type{Gaussian}; m=[0.0], v=mat(1.0)) = Distribution(Multivariate, Gaussian{Moments}, m=m, v=v)

dims(dist::Distribution{V, Gaussian{Moments}}) where V<:VariateType = size(dist.params[:m])

vague(::Type{Gaussian{Moments}}) = Distribution(Univariate, Gaussian{Moments}, m=0.0, v=huge)
vague(::Type{Gaussian{Moments}}, dims::Tuple{Int64}) = Distribution(Multivariate, Gaussian{Moments}, m=zeros(dims), v=huge*diageye(dims[1]))
# Default to Moments parameterization
vague(::Type{Gaussian}) = vague(Gaussian{Moments})
vague(::Type{Gaussian}, dims::Tuple{Int64}) = vague(Gaussian{Moments}, dims)

unsafeMean(dist::Distribution{V, Gaussian{Moments}}) where V<:VariateType = deepcopy(dist.params[:m]) # unsafe mean

unsafeMode(dist::Distribution{V, Gaussian{Moments}}) where V<:VariateType = deepcopy(dist.params[:m]) # unsafe mode

unsafeVar(dist::Distribution{Univariate, Gaussian{Moments}}) = dist.params[:v] # unsafe variance
unsafeVar(dist::Distribution{Multivariate, Gaussian{Moments}}) = diag(dist.params[:v])

unsafeCov(dist::Distribution{V, Gaussian{Moments}}) where V<:VariateType = deepcopy(dist.params[:v]) # unsafe covariance

unsafeMeanCov(dist::Distribution{V, Gaussian{Moments}}) where V<:VariateType = (deepcopy(dist.params[:m]), deepcopy(dist.params[:v]))

unsafeWeightedMean(dist::Distribution{V, Gaussian{Moments}}) where V<:VariateType = cholinv(dist.params[:v])*dist.params[:m]

unsafePrecision(dist::Distribution{V, Gaussian{Moments}}) where V<:VariateType = cholinv(dist.params[:v])

unsafeMeanPrecision(dist::Distribution{V, Gaussian{Moments}}) where V<:VariateType = (deepcopy(dist.params[:m]), cholinv(dist.params[:v]))

logPdf(dist::Distribution{Univariate, Gaussian{Moments}}, x) = -0.5*(log(2pi) + log(dist.params[:v]) + (x-dist.params[:m])^2/dist.params[:v])
logPdf(dist::Distribution{Multivariate, Gaussian{Moments}}, x) = -0.5*(dims(dist)[1]*log(2pi) + logdet(dist.params[:v]) + transpose(x-dist.params[:m])*cholinv(dist.params[:v])*(x-dist.params[:m]))

# Converting from m,v to xi,w would require two separate inversions of the covariance matrix;
# this function ensures only a single inversion is performed
function unsafeWeightedMeanPrecision(dist::Distribution{V, Gaussian{Moments}}) where V<:VariateType
    W = cholinv(dist.params[:v])
    return (W*dist.params[:m], W)
end

isProper(dist::Distribution{Univariate, Gaussian{Moments}}) = (floatmin(Float64) < dist.params[:v] < floatmax(Float64))
isProper(dist::Distribution{Multivariate, Gaussian{Moments}}) = isRoundedPosDef(dist.params[:v])

function ==(t::Distribution{V, Gaussian{Moments}}, u::Distribution{V, Gaussian{Moments}}) where V<:VariateType
    (t === u) && return true
    isApproxEqual(t.params[:m], u.params[:m]) && isApproxEqual(t.params[:v], u.params[:v])
end

function ==(t::Distribution{V, Gaussian{Moments}}, u::Distribution{V, F}) where {V<:VariateType, F<:Gaussian}
    (t === u) && return true
    isApproxEqual(t.params[:m], unsafeMean(u)) && isApproxEqual(t.params[:v], unsafeCov(u))
end

# Average energy functional
function averageEnergy(::Type{Gaussian{Moments}}, marg_out::Distribution{Univariate}, marg_mean::Distribution{Univariate}, marg_var::Distribution{Univariate})
    (m_mean, v_mean) = unsafeMeanCov(marg_mean)
    (m_out, v_out) = unsafeMeanCov(marg_out)

    0.5*log(2*pi) +
    0.5*unsafeLogMean(marg_var) +
    0.5*unsafeInverseMean(marg_var)*(v_out + v_mean + (m_out - m_mean)^2)
end

function averageEnergy(::Type{Gaussian{Moments}}, marg_out::Distribution{Multivariate}, marg_mean::Distribution{Multivariate}, marg_var::Distribution{MatrixVariate})
    (m_mean, v_mean) = unsafeMeanCov(marg_mean)
    (m_out, v_out) = unsafeMeanCov(marg_out)

    0.5*dims(marg_out)[1]*log(2*pi) +
    0.5*unsafeDetLogMean(marg_var) +
    0.5*tr( unsafeInverseMean(marg_var)*(v_out + v_mean + (m_out - m_mean)*(m_out - m_mean)'))
end

function averageEnergy(::Type{Gaussian{Moments}}, marg_out_mean::Distribution{Multivariate, <:Gaussian}, marg_var::Distribution{Univariate})
    (m, V) = unsafeMeanCov(marg_out_mean)

    0.5*log(2*pi) +
    0.5*unsafeLogMean(marg_var) +
    0.5*unsafeInverseMean(marg_var)*( V[1,1] - V[1,2] - V[2,1] + V[2,2] + (m[1] - m[2])^2 )
end

function averageEnergy(::Type{Gaussian{Moments}}, marg_out_mean::Distribution{Multivariate, <:Gaussian}, marg_var::Distribution{MatrixVariate})
    (m, V) = unsafeMeanCov(marg_out_mean)
    d = Int64(dims(marg_out_mean)[1]/2)

    0.5*d*log(2*pi) +
    0.5*unsafeDetLogMean(marg_var) +
    0.5*tr( unsafeInverseMean(marg_var)*( V[1:d,1:d] - V[1:d,d+1:end] - V[d+1:end,1:d] + V[d+1:end,d+1:end] + (m[1:d] - m[d+1:end])*(m[1:d] - m[d+1:end])' ) )
end
