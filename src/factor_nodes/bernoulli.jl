export Bernoulli

"""
Description:
    Bernoulli factor node

    c ∈ {0, 1}
    x ∈ [0, 1]
    
    f(c,x) = Ber(c|x)

Interfaces:
    1. in
    2. out

Construction:
    Bernoulli(id=:some_id)
"""
type Bernoulli <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function Bernoulli(out::Variable, in1::Variable; id=generateId(Bernoulli))
        self = new(id, Array(Interface, 2), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:in] = self.interfaces[1] = associate!(Interface(self), in1)
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