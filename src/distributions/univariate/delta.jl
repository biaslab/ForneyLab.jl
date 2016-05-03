export DeltaDistribution

"""
Description:

    Encodes a delta probability distribution. The DeltaDistribution is used to fix variables to a value, for example to capture observed data.
    p(x) = δ(x-m)

Parameters:

    m (Any)

Construction:

    DeltaDistribution(m)

Example:

    DeltaDistribution(3.0)
"""
type DeltaDistribution{T} <: UnivariateProbabilityDistribution
    m::T
end

DeltaDistribution() = DeltaDistribution{Float64}(1.0)

format(dist::DeltaDistribution) = "δ(m=$(format(dist.m)))"

show(io::IO, dist::DeltaDistribution) = println(io, format(dist))

function prod!{T}(x::DeltaDistribution{T}, y::DeltaDistribution{T}, z::DeltaDistribution{T}=deepcopy(y))
    # Product of two deltas: only valid if the deltas are equal
    (x.m == y.m) || error("The product of two deltas at different positions does not yield a probability distribution")
    (z.m == y.m) || (z.m = deepcopy(y.m))

    return z
end

isProper(dist::DeltaDistribution) = true

Base.mean(dist::DeltaDistribution) = dist.m

Base.var(dist::DeltaDistribution) = 0.0

sample(dist::DeltaDistribution) = dist.m

==(x::DeltaDistribution, y::DeltaDistribution) = (x.m == y.m)

# We can convert a lot of object types into a DeltaDistribution with that object as position of the delta.
# This is useful so we can write i.e. TerminalNode(3.0) instead of TerminalNode(DeltaDistribution(3.0)).
convert(::Type{ProbabilityDistribution}, obj::Number) = DeltaDistribution(obj)

convert(::Type{ProbabilityDistribution}, obj::Bool) = DeltaDistribution(obj)

convert(::Type{DeltaDistribution}, obj::Number) = DeltaDistribution(obj)

convert(::Type{DeltaDistribution}, obj::Bool) = DeltaDistribution(obj)