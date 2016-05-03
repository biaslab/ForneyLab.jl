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

Base.mean(dist::BetaDistribution) = isProper(dist) ? dist.a/(dist.a+dist.b) : NaN

Base.var(dist::BetaDistribution) = isProper(dist) ? dist.a*dist.b/((dist.a+dist.b)^2*(dist.a+dist.b+1)) : NaN

format(dist::BetaDistribution) = "Bet(a=$(format(dist.a)), b=$(format(dist.b)))"

show(io::IO, dist::BetaDistribution) = println(io, format(dist))

function prod!(x::BetaDistribution, y::BetaDistribution, z::BetaDistribution=BetaDistribution())
    # Multiplication of 2 beta PDFs: p(z) = p(x) * p(y)
    z.a = x.a + y.a - 1.0
    z.b = x.b + y.b - 1.0

    return z
end

@symmetrical function prod!(x::BetaDistribution, y::DeltaDistribution{Float64}, z::DeltaDistribution{Float64}=DeltaDistribution(1.0))
    # Multiplication of beta PDF with Delta
    (0.0 <= y.m <= 1.0) || throw(DomainError())
    z.m = y.m

    return z
end

@symmetrical function prod!(x::BetaDistribution, y::DeltaDistribution{Float64}, z::BetaDistribution)
    # Multiplication of beta PDF with Delta, force result to be beta
    (0.0 <= y.m <= 1.0) || throw(DomainError())
    z.a = (1.0 - y.m) * huge
    z.b = y.m * huge

    return z
end

==(x::BetaDistribution, y::BetaDistribution) = (x.a==y.a && x.b==y.b)
