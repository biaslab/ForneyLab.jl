export InverseGammaDistribution

"""
Description:

    Encodes an inverse gamma PDF.

Pamameters:

    Real scalars a > 0 (shape) and b > 0 (scale).

Construction:

    InverseGammaDistribution(a=1.0, b=1.0)

Reference:

    Korl, 2005; A factor graph approach to signal modelling, system identification and filtering; appendix A
"""
type InverseGammaDistribution <: UnivariateProbabilityDistribution
    a::Float64 # shape
    b::Float64 # scale
end

InverseGammaDistribution(; a=3.0, b=2.0) = InverseGammaDistribution(a, b)

function vague!(dist::InverseGammaDistribution)
    dist.a = tiny
    dist.b = tiny
    return dist
end

isProper(dist::InverseGammaDistribution) = (dist.a >= tiny && dist.b >= tiny)

function Base.mean(dist::InverseGammaDistribution)
    if isProper(dist) && dist.a > 1.0
        return dist.b / (dist.a - 1)
    else
        return NaN
    end
end

function Base.var(dist::InverseGammaDistribution)
    if isProper(dist) && dist.a > 2.0
        return (dist.b^2) / ( ( (dist.a-1)^2 ) * (dist.a-2) )
    else
        return NaN
    end
end

format(dist::InverseGammaDistribution) = "Ig(a=$(format(dist.a)), b=$(format(dist.b)))"

show(io::IO, dist::InverseGammaDistribution) = println(io, format(dist))

function prod!(x::InverseGammaDistribution, y::InverseGammaDistribution, z::InverseGammaDistribution=InverseGammaDistribution())
    # Multiplication of 2 inverse gamma PDFs: p(z) = p(x) * p(y)
    z.a = x.a + y.a + 1.0
    z.b = x.b + y.b

    return z
end

@symmetrical function prod!(x::InverseGammaDistribution, y::DeltaDistribution{Float64}, z::DeltaDistribution{Float64}=DeltaDistribution(1.0))
    # Multiplication of inverse gamma PDF with a Delta
    (y.m >= 0.0) || throw(DomainError())
    z.m = y.m

    return z
end

@symmetrical function prod!(x::InverseGammaDistribution, y::DeltaDistribution{Float64}, z::InverseGammaDistribution)
    # Multiplication of inverse gamma PDF with a Delta, force result to be inverse gamma
    (y.m >= 0.0) || throw(DomainError())
    z.a = clamp(1e8*y.m, tiny, huge)
    z.b = clamp(y.m*(z.a-1), tiny, huge)

    return z
end

==(x::InverseGammaDistribution, y::InverseGammaDistribution) = (x.a==y.a && x.b==y.b)