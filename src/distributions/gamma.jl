############################################
# GammaDistribution
############################################
# Description:
#   Encodes a gamma PDF.
#   Pamameters: scalars a (shape) and b (rate).
############################################
export GammaDistribution

type GammaDistribution <: ProbabilityDistribution
    a::Float64 # shape
    b::Float64 # rate
    GammaDistribution(; a=1.0, b=1.0) = new(a, b)
end

uninformative(dist_type::Type{GammaDistribution}) = GammaDistribution(a=0.999, b=0.001)
Base.mean(dist::GammaDistribution) = dist.a / dist.b
Base.var(dist::GammaDistribution) = dist.a / (dist.b^2)

function show(io::IO, dist::GammaDistribution)
    println(io, typeof(dist))
    println(io, "a = $(dist.a) (shape)")
    println(io, "b = $(dist.b) (rate)")
end

function calculateMarginal(forward_dist::GammaDistribution, backward_dist::GammaDistribution)
    return GammaDistribution(a = forward_dist.a+backward_dist.a-1.0, b = forward_dist.b+backward_dist.b)    
end

function calculateMarginal!(edge::Edge, forward_dist::GammaDistribution, backward_dist::GammaDistribution)
    # Calculate the marginal from a forward/backward message pair.
    # We calculate the marginal by using the EqualityNode update rules; same for the functions below
    marg = getOrCreateMarginal(edge, GammaDistribution)
    marg.a = forward_dist.a+backward_dist.a-1.0
    marg.b = forward_dist.b+backward_dist.b
    return marg
end