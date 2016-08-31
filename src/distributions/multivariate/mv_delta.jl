export MvDelta

"""
Description:

    Encodes a multivariate delta distribution. Can be used to carry samples/observations.
    p(x) = δ(x-m)

Parameters:

    m (location) (Any vector)

Construction:

    MvDelta([1.0, 3.0])
"""
type MvDelta{T, dims} <: Multivariate{dims}
    m::Vector{T}
end

MvDelta{T<:Any}(m::Vector{T}) = MvDelta{T, length(m)}(deepcopy(m))

MvDelta() = MvDelta{Float64, 1}(ones(1))

pdf(dist::MvDelta, x::Vector) = ((dist.m==x) ? 1.0 : 0.0)

format(dist::MvDelta) = "δ(m=$(format(dist.m)))"

show(io::IO, dist::MvDelta) = println(io, format(dist))

function prod!{T,dims}(x::MvDelta{T,dims}, y::MvDelta{T,dims}, z::MvDelta{T,dims}=deepcopy(y))
    # Sifting property: multiplying a prob. dist. with a delta yields the delta
    (x.m == y.m) || error("The product of two deltas at different positions does not yield a probability distribution")
    (z.m == y.m) || (z.m = deepcopy(y.m))
    return z
end

@symmetrical function prod!{dims}(::Void, y::MvDelta{Float64,dims}, z::MvDelta{Float64,dims}=MvDelta(y.m))
    # Product of an unknown PDF and MvDelta
    (z.m == y.m) || (z.m[:] = y.m)

    return z
end

isProper(dist::MvDelta) = true

vague{T,dims}(::Type{MvDelta{T,dims}}) = MvDelta(rand(T, dims))

m(dist::MvDelta) = dist.m

V(dist::MvDelta) = zeros(dist.m)

S(dist::MvDelta) = Diagonal(zeros(dist.m))

sample(dist::MvDelta) = dist.m

==(x::MvDelta, y::MvDelta) = (x.m == y.m)

# We can convert a lot of object types into a MvDelta with that object as position of the delta.
# This is useful so we can write i.e. TerminalNode([1.0, 3.0]) instead of TerminalNode(MvDelta([1.0, 3.0])).
convert{T<:Number}(::Type{ProbabilityDistribution}, obj::Vector{T}) = MvDelta(obj)

convert{T<:Number}(::Type{MvDelta}, obj::Vector{T}) = MvDelta(obj)
