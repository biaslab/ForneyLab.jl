export BetaDistribution

"""
Description:

    Encodes a beta PDF.

Pamameters:

    Real scalars a > 0 (shape) and b > 0 (rate).

Construction:

    BetaDistribution(a=1.0, b=1.0)

Reference:

    Bishop, 2006; Pattern recognition and machine learning; appendix B
"""
type BetaDistribution <: UnivariateProbabilityDistribution
    a::Float64 # shape
    b::Float64 # rate
end

BetaDistribution(; a=1.0, b=1.0) = BetaDistribution(a, b)

function vague!(dist::BetaDistribution)
    dist.a = tiny
    dist.b = tiny
    return dist
end

isProper(dist::BetaDistribution) = (dist.a >= tiny && dist.b >= tiny)

Base.mean(dist::BetaDistribution) = isProper(dist) ? a/(a+b) : NaN

Base.var(dist::BetaDistribution) = isProper(dist) ? a*b/((a+b)^2*(a+b+1)) : NaN

format(dist::BetaDistribution) = "Bet(a=$(format(dist.a)), b=$(format(dist.b)))"

show(io::IO, dist::BetaDistribution) = println(io, format(dist))

==(x::BetaDistribution, y::BetaDistribution) = (x.a==y.a && x.b==y.b)
