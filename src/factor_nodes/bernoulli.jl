export Bernoulli

"""
Description:
    Bernoulli factor node

    z ∈ {0, 1}
    p ∈ [0, 1]
    
    f(z,p) = Ber(z|p)

Interfaces:
    1. p
    2. out

Construction:
    Bernoulli(id=:some_id)
"""
type Bernoulli <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function Bernoulli(out::Variable, p::Variable; id=generateId(Bernoulli))
        self = new(id, Array(Interface, 2), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:p] = self.interfaces[1] = associate!(Interface(self), p)
        self.i[:out] = self.interfaces[2] = associate!(Interface(self), out)

        return self
    end
end

slug(::Type{Bernoulli}) = "Ber"

ProbabilityDistribution(::Type{Bernoulli}) = ProbabilityDistribution(Bernoulli, p=0.5)

vague(::Type{ProbabilityDistribution{Bernoulli}}) = ProbabilityDistribution(Bernoulli, p=0.5)

isProper(dist::ProbabilityDistribution{Bernoulli}) = (0 <= dist.params[:p] <= 1)

unsafeMean(dist::ProbabilityDistribution{Bernoulli}) = dist.params[:p]

unsafeVar(dist::ProbabilityDistribution{Bernoulli}) = dist.params[:p]*(1-dist.params[:p])

function prod!( x::ProbabilityDistribution{Bernoulli},
                y::ProbabilityDistribution{Bernoulli},
                z::ProbabilityDistribution{Bernoulli}=ProbabilityDistribution(Bernoulli, p=0.5))

    norm = x.params[:p] * y.params[:p] + (1 - x.params[:p]) * (1 - y.params[:p])
    (norm > 0) || error("Product of $(x) and $(y) cannot be normalized")
    z.params[:p] = (x.params[:p] * y.params[:p]) / norm

    return z
end

# Entropy functional
function differentialEntropy(dist::ProbabilityDistribution{Bernoulli})
    -(1.0 - dist.params[:p])*log(1.0 - dist.params[:p]) -
    dist.params[:p]*log(dist.params[:p])
end

# Average energy functional
function averageEnergy(::Type{Bernoulli}, marg_in::ProbabilityDistribution, marg_out::ProbabilityDistribution)
    -unsafeMean(marg_out)*unsafeLogMean(marg_in) -
    (1.0 - unsafeMean(marg_out))*unsafeMirroredLogMean(marg_in)
end
