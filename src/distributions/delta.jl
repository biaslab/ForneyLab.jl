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
DeltaDistribution() = DeltaDistribution(1.0)

show(io::IO, dist::DeltaDistribution) = println(io, "DeltaDistribution($(dist.m))")

Base.mean(dist::DeltaDistribution) = dist.m
Base.var(dist::DeltaDistribution) = 0.0

==(x::DeltaDistribution, y::DeltaDistribution) = (x.m == y.m)

# We can convert any object into a DeltaDistribution with that object as position of the delta
# This is useful so we can write i.e. TerminalNode(3.0) instead of TerminalNode(DeltaDistribution(3.0))
convert{T<:Any}(::Type{DeltaDistribution{T}}, obj::T) = DeltaDistribution(obj)
convert(::Type{ProbabilityDistribution}, obj::Any) = DeltaDistribution(obj)