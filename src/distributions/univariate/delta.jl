export Delta

"""
Description:

    Encodes a delta probability distribution. The Delta is used to fix variables to a value, for example to capture observed data.
    p(x) = δ(x-m)

Parameters:

    m (Any)

Construction:

    Delta(m)

Example:

    Delta(3.0)
"""
type Delta{T} <: Univariate
    m::T
end

Delta() = Delta{Float64}(1.0)

pdf(dist::Delta, x) = ((dist.m==x) ? 1.0 :0.0)

format(dist::Delta) = "δ(m=$(format(dist.m)))"

show(io::IO, dist::Delta) = println(io, format(dist))

function prod!{T}(x::Delta{T}, y::Delta{T}, z::Delta{T}=deepcopy(y))
    # Product of two deltas: only valid if the deltas are equal
    (x.m == y.m) || error("The product of two deltas at different positions does not yield a probability distribution")
    (z.m == y.m) || (z.m = deepcopy(y.m))

    return z
end

isProper(dist::Delta) = true

Base.mean(dist::Delta) = dist.m

Base.var(dist::Delta) = 0.0

sample(dist::Delta) = dist.m

==(x::Delta, y::Delta) = (x.m == y.m)

vague{T}(::Type{Delta{T}}) = Delta(rand(T))

# We can convert a lot of object types into a Delta with that object as position of the delta.
# This is useful so we can write i.e. TerminalNode(3.0) instead of TerminalNode(Delta(3.0)).
convert(::Type{ProbabilityDistribution}, obj::Number) = Delta(obj)

convert(::Type{ProbabilityDistribution}, obj::Bool) = Delta(obj)

convert(::Type{Delta}, obj::Number) = Delta(obj)

convert(::Type{Delta}, obj::Bool) = Delta(obj)