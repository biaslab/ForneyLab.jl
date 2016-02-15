export BernoulliDistribution

"""
Description:

    Encodes a distribution over binary domain {false,true}.

Pamameters:

    p ∈ [0,1]: Pr{X=true} = p

Construction:

    BernoulliDistribution(p)
"""
type BernoulliDistribution <: UnivariateProbabilityDistribution
    p::Float64 # Pr{X=true}
end
BernoulliDistribution() = BernoulliDistribution(0.5)

function vague!(dist::BernoulliDistribution)
    dist.p = 0.5
    return dist
end

isProper(dist::BernoulliDistribution) = (0 <= dist.p <= 1)

Base.mean(dist::BernoulliDistribution) = dist.p

Base.mean(::Type{DeltaDistribution{Float64}}, d::BernoulliDistribution) = DeltaDistribution(mean(d)) # Definition for post-processing

Base.var(dist::BernoulliDistribution) = dist.p*(1-dist.p)

sample(dist::BernoulliDistribution) = (rand() < dist.p)

sample(::Type{DeltaDistribution{Bool}}, d::BernoulliDistribution) = DeltaDistribution(sample(d)) # Definition for post-processing

format(dist::BernoulliDistribution) = "Bernoulli(p=$(format(dist.p)))"

show(io::IO, dist::BernoulliDistribution) = println(io, format(dist))

==(x::BernoulliDistribution, y::BernoulliDistribution) = isApproxEqual(x.p, y.p)

# Converts from DeltaDistribution{Bool} -> BernoulliDistribution
Base.convert{T<:Bool}(::Type{BernoulliDistribution}, delta::DeltaDistribution{T}) = BernoulliDistribution(float(delta.m))
