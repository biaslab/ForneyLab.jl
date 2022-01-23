format(dist::Distribution{V, GaussianMeanVariance}) where V<:VariateType = "$(slug(GaussianMeanVariance))(m=$(format(dist.params[:m])), v=$(format(dist.params[:v])))"

Distribution(::Type{Univariate}, ::Type{GaussianMeanVariance}; m=0.0, v=1.0) = Distribution{Univariate, GaussianMeanVariance}(Dict(:m=>m, :v=>v))
Distribution(::Type{GaussianMeanVariance}; m::Number=0.0, v::Number=1.0) = Distribution{Univariate, GaussianMeanVariance}(Dict(:m=>m, :v=>v))
Distribution(::Type{Multivariate}, ::Type{GaussianMeanVariance}; m=[0.0], v=mat(1.0)) = Distribution{Multivariate, GaussianMeanVariance}(Dict(:m=>m, :v=>v))

dims(dist::Distribution{V, GaussianMeanVariance}) where V<:VariateType = size(dist.params[:m])

vague(::Type{GaussianMeanVariance}) = Distribution(Univariate, GaussianMeanVariance, m=0.0, v=huge)
vague(::Type{GaussianMeanVariance}, dims::Tuple{Int64}) = Distribution(Multivariate, GaussianMeanVariance, m=zeros(dims), v=huge*diageye(dims[1]))

unsafeMean(dist::Distribution{V, GaussianMeanVariance}) where V<:VariateType = deepcopy(dist.params[:m]) # unsafe mean

unsafeMode(dist::Distribution{V, GaussianMeanVariance}) where V<:VariateType = deepcopy(dist.params[:m]) # unsafe mode

unsafeVar(dist::Distribution{Univariate, GaussianMeanVariance}) = dist.params[:v] # unsafe variance
unsafeVar(dist::Distribution{Multivariate, GaussianMeanVariance}) = diag(dist.params[:v])

unsafeCov(dist::Distribution{V, GaussianMeanVariance}) where V<:VariateType = deepcopy(dist.params[:v]) # unsafe covariance

unsafeMeanCov(dist::Distribution{V, GaussianMeanVariance}) where V<:VariateType = (deepcopy(dist.params[:m]), deepcopy(dist.params[:v]))

unsafeWeightedMean(dist::Distribution{V, GaussianMeanVariance}) where V<:VariateType = cholinv(dist.params[:v])*dist.params[:m]

unsafePrecision(dist::Distribution{V, GaussianMeanVariance}) where V<:VariateType = cholinv(dist.params[:v])

unsafeMeanPrecision(dist::Distribution{V, GaussianMeanVariance}) where V<:VariateType = (deepcopy(dist.params[:m]), cholinv(dist.params[:v]))

logPdf(dist::Distribution{Univariate, GaussianMeanVariance}, x) = -0.5*(log(2pi) + log(dist.params[:v]) + (x-dist.params[:m])^2/dist.params[:v])
logPdf(dist::Distribution{Multivariate, GaussianMeanVariance}, x) = -0.5*(dims(dist)[1]*log(2pi) + logdet(dist.params[:v]) + transpose(x-dist.params[:m])*cholinv(dist.params[:v])*(x-dist.params[:m]))

# Converting from m,v to xi,w would require two separate inversions of the covariance matrix;
# this function ensures only a single inversion is performed
function unsafeWeightedMeanPrecision(dist::Distribution{V, GaussianMeanVariance}) where V<:VariateType
    W = cholinv(dist.params[:v])
    return (W*dist.params[:m], W)
end

isProper(dist::Distribution{Univariate, GaussianMeanVariance}) = (floatmin(Float64) < dist.params[:v] < floatmax(Float64))
isProper(dist::Distribution{Multivariate, GaussianMeanVariance}) = isRoundedPosDef(dist.params[:v])

function ==(t::Distribution{V, GaussianMeanVariance}, u::Distribution{V, GaussianMeanVariance}) where V<:VariateType
    (t === u) && return true
    isApproxEqual(t.params[:m], u.params[:m]) && isApproxEqual(t.params[:v], u.params[:v])
end

function ==(t::Distribution{V, GaussianMeanVariance}, u::Distribution{V, F}) where {V<:VariateType, F<:Gaussian}
    (t === u) && return true
    isApproxEqual(t.params[:m], unsafeMean(u)) && isApproxEqual(t.params[:v], unsafeCov(u))
end

# Average energy functional
function averageEnergy(::Type{GaussianMeanVariance}, marg_out::Distribution{Univariate}, marg_mean::Distribution{Univariate}, marg_var::Distribution{Univariate})
    (m_mean, v_mean) = unsafeMeanCov(marg_mean)
    (m_out, v_out) = unsafeMeanCov(marg_out)

    0.5*log(2*pi) +
    0.5*unsafeLogMean(marg_var) +
    0.5*unsafeInverseMean(marg_var)*(v_out + v_mean + (m_out - m_mean)^2)
end

function averageEnergy(::Type{GaussianMeanVariance}, marg_out::Distribution{Multivariate}, marg_mean::Distribution{Multivariate}, marg_var::Distribution{MatrixVariate})
    (m_mean, v_mean) = unsafeMeanCov(marg_mean)
    (m_out, v_out) = unsafeMeanCov(marg_out)

    0.5*dims(marg_out)[1]*log(2*pi) +
    0.5*unsafeDetLogMean(marg_var) +
    0.5*tr( unsafeInverseMean(marg_var)*(v_out + v_mean + (m_out - m_mean)*(m_out - m_mean)'))
end

function averageEnergy(::Type{GaussianMeanVariance}, marg_out_mean::Distribution{Multivariate, <:Gaussian}, marg_var::Distribution{Univariate})
    (m, V) = unsafeMeanCov(marg_out_mean)

    0.5*log(2*pi) +
    0.5*unsafeLogMean(marg_var) +
    0.5*unsafeInverseMean(marg_var)*( V[1,1] - V[1,2] - V[2,1] + V[2,2] + (m[1] - m[2])^2 )
end

function averageEnergy(::Type{GaussianMeanVariance}, marg_out_mean::Distribution{Multivariate, <:Gaussian}, marg_var::Distribution{MatrixVariate})
    (m, V) = unsafeMeanCov(marg_out_mean)
    d = Int64(dims(marg_out_mean)[1]/2)

    0.5*d*log(2*pi) +
    0.5*unsafeDetLogMean(marg_var) +
    0.5*tr( unsafeInverseMean(marg_var)*( V[1:d,1:d] - V[1:d,d+1:end] - V[d+1:end,1:d] + V[d+1:end,d+1:end] + (m[1:d] - m[d+1:end])*(m[1:d] - m[d+1:end])' ) )
end
