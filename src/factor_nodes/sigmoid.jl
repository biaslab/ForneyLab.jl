export Sigmoid

"""
Description:
    Constrains a continuous, real-valued variable with a binary (boolean) variable.

    f(x,y) = σ(x⋅y)

Interfaces:
    1. real
    2. bin

Construction:
    Sigmoid(id=:some_id)
"""
type Sigmoid <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function Sigmoid(bin::Variable, real::Variable; id=generateId(Sigmoid))
        self = new(id, Array(Interface, 2), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:real] = self.interfaces[1] = associate!(Interface(self), real)
        self.i[:bin] = self.interfaces[2] = associate!(Interface(self), bin)

        return self
    end
end

slug(::Type{Sigmoid}) = "σ"

# Average energy functional
function averageEnergy(::Type{Sigmoid}, marg_real::ProbabilityDistribution{Gaussian}, marg_bin::ProbabilityDistribution{Bernoulli})
    ensureParameters!(marg_real, (:m, :v))
    h = (x -> log(0.5*erf(x) + 0.5 + tiny))   
    
    (1 - marg_bin.params[:p])*gaussianQuadrature(h, m=-marg_real.params[:m], v=marg_real.params[:v]) +
    marg_bin.params[:p]*gaussianQuadrature(h, m=marg_real.params[:m], v=marg_real.params[:v])
end

function averageEnergy(::Type{Sigmoid}, marg_real::ProbabilityDistribution{Gaussian}, marg_bin::ProbabilityDistribution{PointMass})
    p = mapToBernoulliParameterRange(marg_bin.params[:m])

    return averageEnergy(Sigmoid, marg_real, ProbabilityDistribution(Bernoulli, p=p))
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