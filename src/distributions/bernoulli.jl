################################################################
# BernoulliDistribution
################################################################
# Description:
#   Encodes a distribution over binary domain {false,true}.
#   Pamameters: p ∈ [0,1], Pr{X=true} = p
################################################################
export BernoulliDistribution

type BernoulliDistribution <: ProbabilityDistribution
    p::Float64 # Pr{X=true}
end
BernoulliDistribution() = BernoulliDistribution(0.5)

vague(::Type{BernoulliDistribution}) = BernoulliDistribution(0.5)

Base.mean(dist::BernoulliDistribution) = dist.p

Base.var(dist::BernoulliDistribution) = dist.p*(1-dist.p)

sample(dist::BernoulliDistribution) = (rand() < dist.p)

format(dist::BernoulliDistribution) = "Bernoulli(p=$(format(dist.p)))"
show(io::IO, dist::BernoulliDistribution) = println(io, format(dist))

==(x::BernoulliDistribution, y::BernoulliDistribution) = (x.p == y.p)

# Converts from DeltaDistribution{Bool} -> BernoulliDistribution
convert{T<:Bool}(::Type{BernoulliDistribution}, delta::DeltaDistribution{T}) = BernoulliDistribution(float(delta.m))
convert{T<:Bool}(::Type{Message{BernoulliDistribution}}, msg::Message{DeltaDistribution{T}}) = Message(BernoulliDistribution(float(delta.m)))
