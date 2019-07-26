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
        self = new(id, Array{Interface}(undef, 3), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:m] = self.interfaces[2] = associate!(Interface(self), m)
        self.i[:w] = self.interfaces[3] = associate!(Interface(self), w)

        return self
    end
end

slug(::Type{GaussianMeanPrecision}) = "𝒩"

function format(dist::ProbabilityDistribution{V, GaussianMeanPrecision}) where V<:VariateType
    return "$(slug(GaussianMeanPrecision))(m=$(format(dist.params[:m])), w=$(format(dist.params[:w])))"
end

function ProbabilityDistribution(::Type{Univariate}, ::Type{GaussianMeanPrecision}; m=0.0, w=1.0)
    return ProbabilityDistribution{Univariate, GaussianMeanPrecision}(Dict(:m=>m, :w=>w))
end

function ProbabilityDistribution(::Type{GaussianMeanPrecision}; m::Number=0.0, w::Number=1.0)
    return ProbabilityDistribution{Univariate, GaussianMeanPrecision}(Dict(:m=>m, :w=>w))
end

function ProbabilityDistribution(::Type{Multivariate}, ::Type{GaussianMeanPrecision}; m=[0.0], w=Matrix{Float64}(I,1,1))
    return ProbabilityDistribution{Multivariate, GaussianMeanPrecision}(Dict(:m=>m, :w=>PDMat(Matrix(w))))
end

function dims(dist::ProbabilityDistribution{V, GaussianMeanPrecision}) where V<:VariateType
    return length(dist.params[:m])
end

function vague(::Type{GaussianMeanPrecision})
    return ProbabilityDistribution(Univariate, GaussianMeanPrecision, m=0.0, w=tiny)
end

function vague(::Type{GaussianMeanPrecision}, dims::Int64)
    w = PDMat(tiny*Matrix{Float64}(I,dims,dims))
    return ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=zeros(dims), w=w)
end

function unsafeMean(dist::ProbabilityDistribution{V, GaussianMeanPrecision}) where V<:VariateType
    return deepcopy(dist.params[:m])
end

function unsafeVar(dist::ProbabilityDistribution{Univariate, GaussianMeanPrecision})
    return 1.0/dist.params[:w]
end

function unsafeVar(dist::ProbabilityDistribution{Multivariate, GaussianMeanPrecision})
    return diag(inv(dist.params[:w]))
end

function unsafeCov(dist::ProbabilityDistribution{V, GaussianMeanPrecision}) where V<:VariateType
    return inv(dist.params[:w])
end

function unsafeMeanCov(dist::ProbabilityDistribution{V, GaussianMeanPrecision}) where V<:VariateType
    return (deepcopy(dist.params[:m]), inv(dist.params[:w]))
end

function unsafeWeightedMean(dist::ProbabilityDistribution{V, GaussianMeanPrecision}) where V<:VariateType
    return dist.params[:w]*dist.params[:m]
end

function unsafePrecision(dist::ProbabilityDistribution{V, GaussianMeanPrecision}) where V<:VariateType
    return deepcopy(dist.params[:w])
end

function unsafeWeightedMeanPrecision(dist::ProbabilityDistribution{V, GaussianMeanPrecision}) where V<:VariateType
    return (dist.params[:w]*dist.params[:m], deepcopy(dist.params[:w]))
end

function isProper(dist::ProbabilityDistribution{Univariate, GaussianMeanPrecision})
    return (floatmin(Float64) < dist.params[:w] < floatmax(Float64))
end

function isProper(dist::ProbabilityDistribution{Multivariate, GaussianMeanPrecision})
    return isRoundedPosDef(dist.params[:w])
end

function ==(t::ProbabilityDistribution{V, GaussianMeanPrecision},
            u::ProbabilityDistribution{V, GaussianMeanPrecision}) where V<:VariateType
    (t === u) && return true
    return isApproxEqual(t.params[:m], u.params[:m]) && isApproxEqual(t.params[:w], u.params[:w])
end

function ==(t::ProbabilityDistribution{V, GaussianMeanPrecision},
            u::ProbabilityDistribution{V, F}) where {V<:VariateType, F<:Gaussian}
    (t === u) && return true
    return isApproxEqual(t.params[:m], unsafeMean(u)) && isApproxEqual(t.params[:w], unsafePrecision(u))
end

# Average energy functional
function averageEnergy(::Type{GaussianMeanPrecision},
                       marg_out::ProbabilityDistribution{Univariate},
                       marg_mean::ProbabilityDistribution{Univariate},
                       marg_prec::ProbabilityDistribution{Univariate})
    (m_mean, v_mean) = unsafeMeanCov(marg_mean)
    (m_out, v_out) = unsafeMeanCov(marg_out)

    0.5*log(2*pi) -
    0.5*unsafeLogMean(marg_prec) +
    0.5*unsafeMean(marg_prec)*(v_out + v_mean + (m_out - m_mean)^2)
end

function averageEnergy(::Type{GaussianMeanPrecision},
                       marg_out::ProbabilityDistribution{Multivariate},
                       marg_mean::ProbabilityDistribution{Multivariate},
                       marg_prec::ProbabilityDistribution{MatrixVariate})
    (m_mean, v_mean) = unsafeMeanCov(marg_mean)
    (m_out, v_out) = unsafeMeanCov(marg_out)

    0.5*dims(marg_out)*log(2*pi) -
    0.5*unsafeDetLogMean(marg_prec) +
    0.5*tr( unsafeMean(marg_prec)*(v_out + v_mean + (m_out - m_mean)*(m_out - m_mean)' ))
end

function averageEnergy(::Type{GaussianMeanPrecision},
                       marg_out_mean::ProbabilityDistribution{Multivariate, F},
                       marg_prec::ProbabilityDistribution{Univariate}) where F<:Gaussian
    (m, V) = unsafeMeanCov(marg_out_mean)

    0.5*log(2*pi) -
    0.5*unsafeLogMean(marg_prec) +
    0.5*unsafeMean(marg_prec)*( V[1,1] - V[1,2] - V[2,1] + V[2,2] + (m[1] - m[2])^2 )
end

function averageEnergy(::Type{GaussianMeanPrecision},
                       marg_out_mean::ProbabilityDistribution{Multivariate, F},
                       marg_prec::ProbabilityDistribution{MatrixVariate}) where F<:Gaussian
    (m, V) = unsafeMeanCov(marg_out_mean)
    d = Int64(dims(marg_out_mean)/2)

    0.5*d*log(2*pi) -
    0.5*unsafeDetLogMean(marg_prec) +
    0.5*tr( unsafeMean(marg_prec)*( V[1:d,1:d] - V[1:d,d+1:end] - V[d+1:end,1:d] + V[d+1:end,d+1:end] + (m[1:d] - m[d+1:end])*(m[1:d] - m[d+1:end])' ) )
end
