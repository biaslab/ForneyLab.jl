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

type DeltaDistribution{T} <: UnivariateProbabilityDistribution
    m::T
end

DeltaDistribution() = DeltaDistribution{Float64}(1.0)

format(dist::DeltaDistribution) = "δ(m=$(format(dist.m)))"

show(io::IO, dist::DeltaDistribution) = println(io, format(dist))

isProper(dist::DeltaDistribution) = true

Base.mean(dist::DeltaDistribution) = dist.m

Base.mean(::Type{DeltaDistribution{Float64}}, d::DeltaDistribution) = deepcopy(d) # Definition for post-processing

Base.var(dist::DeltaDistribution) = 0.0

sample(dist::DeltaDistribution) = dist.m

sample(::Type{DeltaDistribution{Float64}}, d::DeltaDistribution) = deepcopy(d) # Definition for post-processing

==(x::DeltaDistribution, y::DeltaDistribution) = (x.m == y.m)

# We can convert a lot of object types into a DeltaDistribution with that object as position of the delta.
# This is useful so we can write i.e. TerminalNode(3.0) instead of TerminalNode(DeltaDistribution(3.0)).
convert(::Type{ProbabilityDistribution}, obj::Number) = DeltaDistribution(obj)

convert(::Type{ProbabilityDistribution}, obj::Bool) = DeltaDistribution(obj)

convert(::Type{DeltaDistribution}, obj::Number) = DeltaDistribution(obj)

convert(::Type{DeltaDistribution}, obj::Bool) = DeltaDistribution(obj)