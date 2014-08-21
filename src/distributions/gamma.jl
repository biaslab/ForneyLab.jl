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

uninformative(::Type{GammaDistribution}) = GammaDistribution(a=0.999, b=0.001)
Base.mean(dist::GammaDistribution) = dist.a / dist.b
Base.var(dist::GammaDistribution) = dist.a / (dist.b^2)

function show(io::IO, dist::GammaDistribution)
    println(io, typeof(dist))
    println(io, "a = $(dist.a) (shape)")
    println(io, "b = $(dist.b) (rate)")
end

==(x::GammaDistribution, y::GammaDistribution) = (x.a==y.a && x.b==y.b)