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

vague(::Type{GammaDistribution}) = GammaDistribution(a=tiny, b=tiny) # Scale invariant (Jeffrey's) prior
isProper(dist::GammaDistribution) = (dist.a >= tiny && dist.b >= tiny)
Base.mean(dist::GammaDistribution) = isProper(dist) ? dist.a/dist.b : NaN
Base.var(dist::GammaDistribution) = isProper(dist) ? dist.a / (dist.b^2) : NaN
format(dist::GammaDistribution) = "Gam(a=$(format(dist.a)), b=$(format(dist.b)))"
show(io::IO, dist::GammaDistribution) = println(io, format(dist))

==(x::GammaDistribution, y::GammaDistribution) = (x.a==y.a && x.b==y.b)