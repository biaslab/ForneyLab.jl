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

Base.mean(::Type{DeltaDistribution{Float64}}, d::InverseGammaDistribution) = DeltaDistribution(mean(d)) # Definition for post-processing

function Base.var(dist::InverseGammaDistribution)
    if isProper(dist) && dist.a > 2.0
        return (dist.b^2) / ( ( (dist.a-1)^2 ) * (dist.a-2) )
    else
        return NaN
    end
end

format(dist::InverseGammaDistribution) = "Ig(a=$(format(dist.a)), b=$(format(dist.b)))"

show(io::IO, dist::InverseGammaDistribution) = println(io, format(dist))

==(x::InverseGammaDistribution, y::InverseGammaDistribution) = (x.a==y.a && x.b==y.b)