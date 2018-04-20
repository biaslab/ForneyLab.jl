export GaussianMeanPrecision

"""
Description:

    A Gaussian with mean-precision parameterization:

    f(out,m,w) = 𝒩(out|m,w) = (2π)^{-D/2} |w|^{1/2} exp(-1/2 (out - m)' w (out - m))

Interfaces:

    1. out
    2. m (mean)
    3. w (precision)

Construction:

    GaussianMeanPrecision(out, m, w, id=:some_id)
"""
mutable struct GaussianMeanPrecision <: Gaussian
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function GaussianMeanPrecision(out, m, w; id=generateId(GaussianMeanPrecision))
        @ensureVariables(out, m, w)
        self = new(id, Array{Interface}(3), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:m] = self.interfaces[2] = associate!(Interface(self), m)
        self.i[:w] = self.interfaces[3] = associate!(Interface(self), w)

        return self
    end
end

slug(::Type{GaussianMeanPrecision}) = "𝒩"

format{V<:VariateType}(dist::ProbabilityDistribution{V, GaussianMeanPrecision}) = "$(slug(GaussianMeanPrecision))(m=$(format(dist.params[:m])), w=$(format(dist.params[:w])))"

ProbabilityDistribution(::Type{Univariate}, ::Type{GaussianMeanPrecision}; m=0.0, w=1.0) = ProbabilityDistribution{Univariate, GaussianMeanPrecision}(Dict(m=m, w=w))
ProbabilityDistribution(::Type{GaussianMeanPrecision}; m::Number=0.0, w::Number=1.0) = ProbabilityDistribution{Univariate, GaussianMeanPrecision}(Dict(m=m, w=w))
ProbabilityDistribution(::Type{Multivariate}, ::Type{GaussianMeanPrecision}; m=[0.0], w=[1.0].') = ProbabilityDistribution{Multivariate, GaussianMeanPrecision}(Dict(m=m, w=w))

dims{V<:VariateType}(dist::ProbabilityDistribution{V, GaussianMeanPrecision}) = length(dist.params[:m])

vague(::Type{GaussianMeanPrecision}) = ProbabilityDistribution(Univariate, GaussianMeanPrecision, m=0.0, w=tiny)
vague(::Type{GaussianMeanPrecision}, dims::Int64) = ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=zeros(dims), w=tiny*diageye(dims))

unsafeMean{V<:VariateType}(dist::ProbabilityDistribution{V, GaussianMeanPrecision}) = deepcopy(dist.params[:m]) # unsafe mean

unsafeVar(dist::ProbabilityDistribution{Univariate, GaussianMeanPrecision}) = 1.0/dist.params[:w] # unsafe variance
unsafeVar(dist::ProbabilityDistribution{Multivariate, GaussianMeanPrecision}) = diag(cholinv(dist.params[:w]))

unsafeCov(dist::ProbabilityDistribution{Univariate, GaussianMeanPrecision}) = 1.0/dist.params[:w] # unsafe covariance
unsafeCov(dist::ProbabilityDistribution{Multivariate, GaussianMeanPrecision}) = cholinv(dist.params[:w])

unsafeWeightedMean{V<:VariateType}(dist::ProbabilityDistribution{V, GaussianMeanPrecision}) = dist.params[:w]*dist.params[:m] # unsafe weighted mean

unsafePrecision{V<:VariateType}(dist::ProbabilityDistribution{V, GaussianMeanPrecision}) = deepcopy(dist.params[:w]) # unsafe precision

isProper(dist::ProbabilityDistribution{Univariate, GaussianMeanPrecision}) = (realmin(Float64) < dist.params[:w] < realmax(Float64))
isProper(dist::ProbabilityDistribution{Multivariate, GaussianMeanPrecision}) = isRoundedPosDef(dist.params[:w])

function =={V<:VariateType}(t::ProbabilityDistribution{V, GaussianMeanPrecision}, u::ProbabilityDistribution{V, GaussianMeanPrecision})
    (t === u) && return true
    isApproxEqual(t.params[:m], u.params[:m]) && isApproxEqual(t.params[:w], u.params[:w])
end

function =={V<:VariateType, F<:Gaussian}(t::ProbabilityDistribution{V, GaussianMeanPrecision}, u::ProbabilityDistribution{V, F})
    (t === u) && return true
    isApproxEqual(t.params[:m], unsafeMean(u)) && isApproxEqual(t.params[:w], unsafePrecision(u))
end

# Average energy functional
function averageEnergy(::Type{GaussianMeanPrecision}, marg_out::ProbabilityDistribution{Univariate}, marg_mean::ProbabilityDistribution{Univariate}, marg_prec::ProbabilityDistribution{Univariate})
    0.5*log(2*pi) -
    0.5*unsafeLogMean(marg_prec) +
    0.5*unsafeMean(marg_prec)*( unsafeCov(marg_out) + unsafeCov(marg_mean) + (unsafeMean(marg_out) - unsafeMean(marg_mean))^2 )
end

function averageEnergy(::Type{GaussianMeanPrecision}, marg_out::ProbabilityDistribution{Multivariate}, marg_mean::ProbabilityDistribution{Multivariate}, marg_prec::ProbabilityDistribution{MatrixVariate})
    0.5*dims(marg_out)*log(2*pi) -
    0.5*unsafeDetLogMean(marg_prec) +
    0.5*trace( unsafeMean(marg_prec)*(unsafeCov(marg_out) + unsafeCov(marg_mean) + (unsafeMean(marg_out) - unsafeMean(marg_mean))*(unsafeMean(marg_out) - unsafeMean(marg_mean))' ))
end

function averageEnergy{F<:Gaussian}(::Type{GaussianMeanPrecision}, marg_out_mean::ProbabilityDistribution{Multivariate, F}, marg_prec::ProbabilityDistribution{Univariate})
    V = unsafeCov(marg_out_mean)
    m = unsafeMean(marg_out_mean)

    0.5*log(2*pi) -
    0.5*unsafeLogMean(marg_prec) +
    0.5*unsafeMean(marg_prec)*( V[1,1] - V[1,2] - V[2,1] + V[2,2] + (m[1] - m[2])^2 )
end

function averageEnergy{F<:Gaussian}(::Type{GaussianMeanPrecision}, marg_out_mean::ProbabilityDistribution{Multivariate, F}, marg_prec::ProbabilityDistribution{MatrixVariate})
    V = unsafeCov(marg_out_mean)
    m = unsafeMean(marg_out_mean)
    d = Int64(dims(marg_out_mean)/2)

    0.5*d*log(2*pi) -
    0.5*unsafeDetLogMean(marg_prec) +
    0.5*trace( unsafeMean(marg_prec)*( V[1:d,1:d] - V[1:d,d+1:end] - V[d+1:end,1:d] + V[d+1:end,d+1:end] + (m[1:d] - m[d+1:end])*(m[1:d] - m[d+1:end])' ) )
end