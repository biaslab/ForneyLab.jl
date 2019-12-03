export GaussianMeanVariance

abstract type Gaussian <: SoftFactor end

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
        self = new(id, Array{Interface}(undef, 3), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:m] = self.interfaces[2] = associate!(Interface(self), m)
        self.i[:v] = self.interfaces[3] = associate!(Interface(self), v)

        return self
    end
end

slug(::Type{GaussianMeanVariance}) = "𝒩"

function format(dist::ProbabilityDistribution{V, GaussianMeanVariance}) where V<:VariateType
    return "$(slug(GaussianMeanVariance))(m=$(format(dist.params[:m])), v=$(format(dist.params[:v])))"
end

function ProbabilityDistribution(::Type{Univariate}, ::Type{GaussianMeanVariance}; m=0.0, v=1.0)
    return ProbabilityDistribution{Univariate, GaussianMeanVariance}(Dict(:m=>m, :v=>v))
end

function ProbabilityDistribution(::Type{GaussianMeanVariance}; m::Number=0.0, v::Number=1.0)
    return ProbabilityDistribution{Univariate, GaussianMeanVariance}(Dict(:m=>m, :v=>v))
end

function ProbabilityDistribution(::Type{Multivariate}, ::Type{GaussianMeanVariance}; m=[0.0], v=Matrix{Float64}(I,1,1))
    if isa(v, PDMat)
        return ProbabilityDistribution{Multivariate, GaussianMeanVariance}(Dict(:m=>m, :v=>v))
    else
        return ProbabilityDistribution{Multivariate, GaussianMeanVariance}(Dict(:m=>m, :v=>PDMat(Matrix(v))))
    end
end

function logPdf(dist::ProbabilityDistribution{Univariate, GaussianMeanVariance},x)
    return -0.5*(log(2pi)+log(dist.params[:v]) + (x-dist.params[:m])^2/dist.params[:v])
end

function logPdf(dist::ProbabilityDistribution{Multivariate, GaussianMeanVariance},x)
    return -0.5*(dims(dist)*log(2pi) + logdet(dist.params[:v]) + transpose(x-dist.params[:m])*inv(dist.params[:v])*(x-dist.params[:m]))
end

function dims(dist::ProbabilityDistribution{V, GaussianMeanVariance}) where V<:VariateType
    return length(dist.params[:m])
end

function vague(::Type{GaussianMeanVariance})
    return ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=huge)
end

function vague(::Type{GaussianMeanVariance}, dims::Int64)
    return ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=zeros(dims), v=PDMat(huge*Matrix{Float64}(I,dims,dims)))
end

function vague(::Type{GaussianMeanVariance}, dims::Tuple{Int64})
    return ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=zeros(dims), v=PDMat(huge*Matrix{Float64}(I,dims[1],dims[1])))
end

function unsafeMode(dist::ProbabilityDistribution{V, GaussianMeanVariance}) where V<:VariateType
    return deepcopy(dist.params[:m])
end

function unsafeMean(dist::ProbabilityDistribution{V, GaussianMeanVariance}) where V<:VariateType
    return deepcopy(dist.params[:m])
end

function unsafeVar(dist::ProbabilityDistribution{Univariate, GaussianMeanVariance})
    return dist.params[:v]
end

function unsafeVar(dist::ProbabilityDistribution{Multivariate, GaussianMeanVariance})
    return diag(dist.params[:v])
end

function unsafePrecision(dist::ProbabilityDistribution{Univariate, GaussianMeanVariance})
    return inv(dist.params[:v])
end

function unsafePrecision(dist::ProbabilityDistribution{Multivariate, GaussianMeanVariance})
    return inv(dist.params[:v]).mat
end

function unsafeCov(dist::ProbabilityDistribution{V, GaussianMeanVariance}) where V<:VariateType
    return deepcopy(dist.params[:v])
end

function unsafeMeanCov(dist::ProbabilityDistribution{V, GaussianMeanVariance}) where V<:VariateType
    return (deepcopy(dist.params[:m]), deepcopy(dist.params[:v]))
end

function unsafeWeightedMean(dist::ProbabilityDistribution{V, GaussianMeanVariance}) where V<:VariateType
    return inv(dist.params[:v])*dist.params[:m]
end

function unsafeWeightedMeanPrecision(dist::ProbabilityDistribution{V, GaussianMeanVariance}) where V<:VariateType
    # Converting from m,v to xi,w would require two separate inversions of the covariance matrix; this function ensures only a singly inversion is performed
    W = inv(dist.params[:v])
    return (W*dist.params[:m], W)
end

function isProper(dist::ProbabilityDistribution{Univariate, GaussianMeanVariance})
    return (floatmin(Float64) < dist.params[:v] < floatmax(Float64))
end

function isProper(dist::ProbabilityDistribution{Multivariate, GaussianMeanVariance})
    return isRoundedPosDef(dist.params[:v])
end

function ==(t::ProbabilityDistribution{V, GaussianMeanVariance},
            u::ProbabilityDistribution{V, GaussianMeanVariance}) where V<:VariateType
    (t === u) && return true
    isApproxEqual(t.params[:m], u.params[:m]) && isApproxEqual(t.params[:v], u.params[:v])
end

function ==(t::ProbabilityDistribution{V, GaussianMeanVariance},
            u::ProbabilityDistribution{V, F}) where {V<:VariateType, F<:Gaussian}
    (t === u) && return true
    isApproxEqual(t.params[:m], unsafeMean(u)) && isApproxEqual(t.params[:v], unsafeCov(u))
end

# Average energy functional
function averageEnergy(::Type{GaussianMeanVariance},
                       marg_out::ProbabilityDistribution{Univariate},
                       marg_mean::ProbabilityDistribution{Univariate},
                       marg_var::ProbabilityDistribution{Univariate})

    (m_mean, v_mean) = unsafeMeanCov(marg_mean)
    (m_out, v_out) = unsafeMeanCov(marg_out)

    0.5*log(2*pi) +
    0.5*unsafeLogMean(marg_var) +
    0.5*unsafeInverseMean(marg_var)*(v_out + v_mean + (m_out - m_mean)^2)
end

function averageEnergy(::Type{GaussianMeanVariance},
                       marg_out::ProbabilityDistribution{Multivariate},
                       marg_mean::ProbabilityDistribution{Multivariate},
                       marg_var::ProbabilityDistribution{MatrixVariate})
    (m_mean, v_mean) = unsafeMeanCov(marg_mean)
    (m_out, v_out) = unsafeMeanCov(marg_out)

    0.5*dims(marg_out)*log(2*pi) +
    0.5*unsafeDetLogMean(marg_var) +
    0.5*tr( unsafeInverseMean(marg_var)*(v_out + v_mean + (m_out - m_mean)*(m_out - m_mean)'))
end
