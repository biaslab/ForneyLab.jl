export GammaDistribution

"""
Description:

    Encodes a gamma PDF.

Pamameters:

	Real scalars scalars a > 0 (shape) and b > 0 (rate).

Construction:

	GammaDistribution(a=1.0, b=1.0)

Reference:

	Bishop, 2006; Pattern recognition and machine learning; appendix B
"""
type GammaDistribution <: UnivariateProbabilityDistribution
    a::Float64 # shape
    b::Float64 # rate
end

GammaDistribution(; a=1.0, b=1.0) = GammaDistribution(a, b)

function vague!(dist::GammaDistribution)
    dist.a = tiny
    dist.b = tiny
    return dist
end

isProper(dist::GammaDistribution) = (dist.a >= tiny && dist.b >= tiny)

Base.mean(dist::GammaDistribution) = isProper(dist) ? dist.a/dist.b : NaN

Base.mean(::Type{DeltaDistribution{Float64}}, d::GammaDistribution) = DeltaDistribution(mean(d)) # Definition for post-processing

Base.var(dist::GammaDistribution) = isProper(dist) ? dist.a / (dist.b^2) : NaN

format(dist::GammaDistribution) = "Gam(a=$(format(dist.a)), b=$(format(dist.b)))"

show(io::IO, dist::GammaDistribution) = println(io, format(dist))

==(x::GammaDistribution, y::GammaDistribution) = (x.a==y.a && x.b==y.b)