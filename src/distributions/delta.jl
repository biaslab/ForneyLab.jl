############################################
# DeltaDistribution
############################################
# Description:
#   Encodes a delta probability distribution.
#   p(x) = 1 if x==m
#        = 0 otherwise
#   Can be used to carry samples/observations.
#   Example: DeltaDistribution(3.0)
############################################

export DeltaDistribution

type DeltaDistribution{T} <: ProbabilityDistribution
    m::T

    DeltaDistribution{T}(m::T) = new(m)
    DeltaDistribution() = new()
end
DeltaDistribution{T}(m::T) = DeltaDistribution{T}(m)

show(io::IO, dist::DeltaDistribution) = println(io, "DeltaDistribution($(dist.m))")

Base.mean(dist::DeltaDistribution) = dist.m
Base.var(dist::DeltaDistribution) = 0.0

==(x::DeltaDistribution, y::DeltaDistribution) = (x.m == y.m)