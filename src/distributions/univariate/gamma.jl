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

Base.var(dist::GammaDistribution) = isProper(dist) ? dist.a / (dist.b^2) : NaN

format(dist::GammaDistribution) = "Gam(a=$(format(dist.a)), b=$(format(dist.b)))"

show(io::IO, dist::GammaDistribution) = println(io, format(dist))

function prod!(x::GammaDistribution, y::GammaDistribution, z::GammaDistribution=GammaDistribution())
    # Multiplication of 2 gamma PDFs: p(z) = p(x) * p(y)
    z.a = x.a + y.a - 1.0
    z.b = x.b + y.b

    return z
end

@symmetrical function prod!(x::GammaDistribution, y::DeltaDistribution{Float64}, z::DeltaDistribution{Float64}=DeltaDistribution(1.0))
    # Multiplication of Gamma PDF with a Delta
    (y.m >= 0.0) || throw(DomainError())
    z.m = y.m

    return z
end

@symmetrical function prod!(x::GammaDistribution, y::DeltaDistribution{Float64}, z::GammaDistribution)
    # Multiplication of Gamma PDF with a Delta, force result to be Gamma
    (y.m >= 0.0) || throw(DomainError())
    z.b = clamp(y.m / 1e-10, tiny, huge)
    z.a = clamp(y.m * z.b, tiny, huge)

    return z
end

==(x::GammaDistribution, y::GammaDistribution) = (x.a==y.a && x.b==y.b)