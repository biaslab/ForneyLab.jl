export MvDeltaDistribution

"""
Description:

    Encodes a multivariate delta distribution. Can be used to carry samples/observations.
    p(x) = δ(x-m)

Parameters:

	m (location) (Any vector)
    
Construction:
    
    MvDeltaDistribution([1.0, 3.0])
"""
type MvDeltaDistribution{T, dims} <: MultivariateProbabilityDistribution
    m::Vector{T}
end

MvDeltaDistribution{T<:Any}(m::Vector{T}) = MvDeltaDistribution{T, length(m)}(deepcopy(m))

MvDeltaDistribution() = MvDeltaDistribution{Float64, 1}(ones(1))

format(dist::MvDeltaDistribution) = "δ(m=$(format(dist.m)))"

show(io::IO, dist::MvDeltaDistribution) = println(io, format(dist))

isProper(dist::MvDeltaDistribution) = true

Base.mean(dist::MvDeltaDistribution) = dist.m

Base.mean(::Type{MvDeltaDistribution{Float64}}, d::MvDeltaDistribution) = deepcopy(d) # Definition for post-processing

Base.var(dist::MvDeltaDistribution) = zeros(length(dist.m))

Base.cov(dist::MvDeltaDistribution) = zeros(length(dist.m), length(dist.m))

sample(dist::MvDeltaDistribution) = dist.m

sample(::Type{MvDeltaDistribution{Float64}}, d::MvDeltaDistribution) = deepcopy(d) # Definition for post-processing

==(x::MvDeltaDistribution, y::MvDeltaDistribution) = (x.m == y.m)

# We can convert a lot of object types into a MvDeltaDistribution with that object as position of the delta.
# This is useful so we can write i.e. TerminalNode([1.0, 3.0]) instead of TerminalNode(MvDeltaDistribution([1.0, 3.0])).
convert{T<:Number}(::Type{ProbabilityDistribution}, obj::Vector{T}) = MvDeltaDistribution(obj)

convert{T<:Number}(::Type{MvDeltaDistribution}, obj::Vector{T}) = MvDeltaDistribution(obj)
