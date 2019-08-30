export Sigmoid

"""
Description:
    Constrains a continuous, real-valued variable with a binary (boolean) variable.

    f(bin, real) = σ(bin⋅real)

Interfaces:

    1. bin
    2. real

Construction:

    Sigmoid(id=:some_id)
"""
mutable struct Sigmoid <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function Sigmoid(bin, real; id=generateId(Sigmoid))
        @ensureVariables(bin, real)
        self = new(id, Array{Interface}(undef, 2), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:bin] = self.interfaces[1] = associate!(Interface(self), bin)
        self.i[:real] = self.interfaces[2] = associate!(Interface(self), real)

        return self
    end
end

slug(::Type{Sigmoid}) = "σ"

# Average energy functional
function averageEnergy(::Type{Sigmoid}, marg_bin::ProbabilityDistribution{Univariate, Bernoulli}, marg_real::ProbabilityDistribution{Univariate, F}) where F<:Gaussian
    (marg_real_m, marg_real_v) = unsafeMeanCov(marg_real)
    h = (x -> log(0.5*erf(x) + 0.5 + tiny)) # Add `tiny` for numeric stability

    (1 - marg_bin.params[:p])*gaussianQuadrature(h, m=-marg_real_m, v=marg_real_v) +
    marg_bin.params[:p]*gaussianQuadrature(h, m=marg_real_m, v=marg_real_v)
end

function averageEnergy(::Type{Sigmoid}, marg_bin::ProbabilityDistribution{Univariate, PointMass}, marg_real::ProbabilityDistribution{Univariate, F}) where F<:Gaussian
    p = mapToBernoulliParameterRange(marg_bin.params[:m])

    return averageEnergy(Sigmoid, ProbabilityDistribution(Univariate, Bernoulli, p=p), marg_real)
end

"""
Map `m` to range of Bernoulli parameter `p` ∈ [0, 1]
"""
function mapToBernoulliParameterRange(m)
    if isa(m, Bool) && (m == true)
        p = 1.0
    elseif isa(m, Bool) && (m == false)
        p = 0.0
    elseif isnan(m)
        p = 0.5
    else
        (-1.0 <= m <= 1.0) || error("Value $m can not be converted to range [0, 1]")
        p = 0.5*m + 0.5
    end

    return p
end