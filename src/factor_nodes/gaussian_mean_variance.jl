export GaussianMeanVariance

"""
Description:

    A Gaussian with mean-variance parameterization:

    f(out,m,v) = 𝒩(out|m,v) = (2π)^{-D/2} |v|^{-1/2} exp(-1/2 (out - m)' v^{-1} (out - m))

Interfaces:

    1. out
    2. m (mean)
    3. v (covariance)

Construction:

    GaussianMeanVariance(out, m, v, id=:some_id)
"""
mutable struct GaussianMeanVariance <: Gaussian
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function GaussianMeanVariance(out, m, v; id=generateId(GaussianMeanVariance))
        @ensureVariables(out, m, v)
        self = new(id, Array{Interface}(3), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:m] = self.interfaces[2] = associate!(Interface(self), m)
        self.i[:v] = self.interfaces[3] = associate!(Interface(self), v)

        return self
    end
end

slug(::Type{GaussianMeanVariance}) = "𝒩"

format{V<:VariateType}(dist::ProbabilityDistribution{V, GaussianMeanVariance}) = "$(slug(GaussianMeanVariance))(m=$(format(dist.params[:m])), v=$(format(dist.params[:v])))"

ProbabilityDistribution(::Type{Univariate}, ::Type{GaussianMeanVariance}; m=0.0, v=1.0) = ProbabilityDistribution{Univariate, GaussianMeanVariance}(Dict(:m=>m, :v=>v))
ProbabilityDistribution(::Type{GaussianMeanVariance}; m::Number=0.0, v::Number=1.0) = ProbabilityDistribution{Univariate, GaussianMeanVariance}(Dict(:m=>m, :v=>v))
ProbabilityDistribution(::Type{Multivariate}, ::Type{GaussianMeanVariance}; m=[0.0], v=mat(1.0)) = ProbabilityDistribution{Multivariate, GaussianMeanVariance}(Dict(:m=>m, :v=>v))

dims{V<:VariateType}(dist::ProbabilityDistribution{V, GaussianMeanVariance}) = length(dist.params[:m])

vague(::Type{GaussianMeanVariance}) = ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=huge)
vague(::Type{GaussianMeanVariance}, dims::Int64) = ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=zeros(dims), v=huge*diageye(dims))

unsafeMean{V<:VariateType}(dist::ProbabilityDistribution{V, GaussianMeanVariance}) = deepcopy(dist.params[:m]) # unsafe mean

unsafeVar(dist::ProbabilityDistribution{Univariate, GaussianMeanVariance}) = dist.params[:v] # unsafe variance
unsafeVar(dist::ProbabilityDistribution{Multivariate, GaussianMeanVariance}) = diag(dist.params[:v])

unsafeCov{V<:VariateType}(dist::ProbabilityDistribution{V, GaussianMeanVariance}) = deepcopy(dist.params[:v]) # unsafe covariance

unsafeWeightedMean(dist::ProbabilityDistribution{Univariate, GaussianMeanVariance}) = dist.params[:m]/dist.params[:v] # unsafe weighted mean
unsafeWeightedMean(dist::ProbabilityDistribution{Multivariate, GaussianMeanVariance}) = cholinv(dist.params[:v])*dist.params[:m]

unsafePrecision(dist::ProbabilityDistribution{Univariate, GaussianMeanVariance}) = 1/dist.params[:v] # unsafe precision
unsafePrecision(dist::ProbabilityDistribution{Multivariate, GaussianMeanVariance}) = cholinv(dist.params[:v])

isProper(dist::ProbabilityDistribution{Univariate, GaussianMeanVariance}) = (realmin(Float64) < dist.params[:v] < realmax(Float64))
isProper(dist::ProbabilityDistribution{Multivariate, GaussianMeanVariance}) = isRoundedPosDef(dist.params[:v])

function =={V<:VariateType}(t::ProbabilityDistribution{V, GaussianMeanVariance}, u::ProbabilityDistribution{V, GaussianMeanVariance})
    (t === u) && return true
    isApproxEqual(t.params[:m], u.params[:m]) && isApproxEqual(t.params[:v], u.params[:v])
end

function =={V<:VariateType, F<:Gaussian}(t::ProbabilityDistribution{V, GaussianMeanVariance}, u::ProbabilityDistribution{V, F})
    (t === u) && return true
    isApproxEqual(t.params[:m], unsafeMean(u)) && isApproxEqual(t.params[:v], unsafeCov(u))
end

# Average energy functional
function averageEnergy(::Type{GaussianMeanVariance}, marg_out::ProbabilityDistribution{Univariate}, marg_mean::ProbabilityDistribution{Univariate}, marg_var::ProbabilityDistribution{Univariate})
    0.5*log(2*pi) +
    0.5*unsafeLogMean(marg_var) +
    0.5*unsafeInverseMean(marg_var)*( unsafeCov(marg_out) + unsafeCov(marg_mean) + (unsafeMean(marg_out) - unsafeMean(marg_mean))^2 )
end

function averageEnergy(::Type{GaussianMeanVariance}, marg_out::ProbabilityDistribution{Multivariate}, marg_mean::ProbabilityDistribution{Multivariate}, marg_var::ProbabilityDistribution{MatrixVariate})
    0.5*dims(marg_out)*log(2*pi) +
    0.5*unsafeDetLogMean(marg_var) +
    0.5*trace( unsafeInverseMean(marg_var)*(unsafeCov(marg_out) + unsafeCov(marg_mean) + (unsafeMean(marg_out) - unsafeMean(marg_mean))*(unsafeMean(marg_out) - unsafeMean(marg_mean))' ))
end