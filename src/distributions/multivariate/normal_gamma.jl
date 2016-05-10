export NormalGammaDistribution

"""
Description:

    Encodes a normal-gamma distribution (bivariate).

Pamameters:

    Real scalars m (location), beta (precision) scalars a (shape) and b (rate).

Construction:

    MvNormalGammaDistribution(m=0.0, beta=1.0, a=1.0, b=1.0)

Reference:

    Bishop, 2006; Pattern recognition and machine learning; appendix B
"""
type NormalGammaDistribution <: MultivariateProbabilityDistribution
    # All univariate, so parameters are floats
    m::Float64    # location
    beta::Float64 # precision
    a::Float64    # shape
    b::Float64    # rate
end

NormalGammaDistribution(; m=0.0, beta=1.0, a=1.0, b=1.0) = NormalGammaDistribution(m, beta, a, b)

# TODO: reference
function vague!(dist::NormalGammaDistribution)
    dist.m = 0.0
    dist.beta = 1.0
    dist.a = tiny
    dist.b = tiny
    return dist
end

vague(::Type{NormalGammaDistribution}) = NormalGammaDistribution(m=0.0, beta=1.0, a=tiny, b=tiny)

isProper(dist::NormalGammaDistribution) = ((dist.beta >= tiny) && (dist.a >= tiny)  && (dist.b >= tiny))

format(dist::NormalGammaDistribution) = "Ng(m=$(format(dist.m)), β=$(format(dist.beta)), a=$(format(dist.a)), b=$(format(dist.b)))"

show(io::IO, dist::NormalGammaDistribution) = println(io, format(dist))

==(x::NormalGammaDistribution, y::NormalGammaDistribution) = (x.m==y.m && x.beta==y.beta && x.a==y.a && x.b==y.b)

Base.mean(dist::NormalGammaDistribution) = [dist.m; dist.a/dist.b]

Base.var(dist::NormalGammaDistribution) = [dist.b/(dist.beta*(dist.a-1)); dist.a/(dist.b^2)]

dimensions(distribution::NormalGammaDistribution) = 2

dimensions(distribution_type::Type{NormalGammaDistribution}) = 2

@symmetrical function prod!(x::NormalGammaDistribution, y::MvDeltaDistribution{Float64,2}, z::MvDeltaDistribution{Float64,2}=MvDeltaDistribution(y.m))
    # Product of normal-gamma PDF and MvDelta
    (y.m[2] >= 0.) || throw(DomainError())
    (z.m == y.m) || (z.m[:] = y.m)

    return z
end

@symmetrical function prod!(x::NormalGammaDistribution, y::MvDeltaDistribution{Float64,2}, z::NormalGammaDistribution)
    # Product of normal-gamma PDF and MvDelta, force result to be mv normal-gamma
    (y.m[2] >= 0.) || throw(DomainError())
    z.b = clamp(y.m[2] / 1e-10, tiny, huge)
    z.a = clamp(y.m[2] * z.b, tiny, huge)
    z.m = y.m[1]
    z.beta = huge

    return z
end
