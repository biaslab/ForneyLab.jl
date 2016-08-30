export NormalGamma

"""
Description:

    Encodes a normal-gamma distribution (bivariate).

Pamameters:

    Real scalars m (location), beta (precision) scalars a (shape) and b (rate).

Construction:

    MvNormalGamma(m=0.0, beta=1.0, a=1.0, b=1.0)

Reference:

    Bishop, 2006; Pattern recognition and machine learning; appendix B
"""
type NormalGamma <: Multivariate
    # All univariate, so parameters are floats
    m::Float64    # location
    beta::Float64 # precision
    a::Float64    # shape
    b::Float64    # rate
end

NormalGamma(; m=0.0, beta=1.0, a=1.0, b=1.0) = NormalGamma(m, beta, a, b)

function pdf(dist::NormalGamma, x::Vector{Float64})
    (length(x) == 2) || return 0.0
    τ = x[2] # the gamma-distributed component
    (τ >= 0) || return 0.0
    m = dist.m; λ = dist.beta; a = dist.a; b = dist.b
    C = b^a * sqrt(λ) / (gamma(a)*sqrt(2.0*pi))

    return (C * τ^(a-0.5) * exp(-b*τ) * exp(-0.5*λ*τ*(x[1]-m)^2))[1]
end

# TODO: reference
function vague!(dist::NormalGamma)
    dist.m = 0.0
    dist.beta = 1.0
    dist.a = tiny
    dist.b = tiny
    return dist
end

vague(::Type{NormalGamma}) = NormalGamma(m=0.0, beta=1.0, a=tiny, b=tiny)

isProper(dist::NormalGamma) = ((dist.beta >= tiny) && (dist.a >= tiny)  && (dist.b >= tiny))

format(dist::NormalGamma) = "Ng(m=$(format(dist.m)), β=$(format(dist.beta)), a=$(format(dist.a)), b=$(format(dist.b)))"

show(io::IO, dist::NormalGamma) = println(io, format(dist))

==(x::NormalGamma, y::NormalGamma) = (x.m==y.m && x.beta==y.beta && x.a==y.a && x.b==y.b)

unsafeMean(dist::NormalGamma) = [dist.m; dist.a/dist.b]

unsafeVar(dist::NormalGamma) = [dist.b/(dist.beta*(dist.a-1)); dist.a/(dist.b^2)]

dimensions(distribution::NormalGamma) = 2

dimensions(distribution_type::Type{NormalGamma}) = 2

@symmetrical function prod!(x::NormalGamma, y::MvDelta{Float64,2}, z::MvDelta{Float64,2}=MvDelta(y.m))
    # Product of normal-gamma PDF and MvDelta
    (y.m[2] >= 0.) || throw(DomainError())
    (z.m == y.m) || (z.m[:] = y.m)

    return z
end

@symmetrical function prod!(x::NormalGamma, y::MvDelta{Float64,2}, z::NormalGamma)
    # Product of normal-gamma PDF and MvDelta, force result to be mv normal-gamma
    (y.m[2] >= 0.) || throw(DomainError())
    z.b = clamp(y.m[2] / 1e-10, tiny, huge)
    z.a = clamp(y.m[2] * z.b, tiny, huge)
    z.m = y.m[1]
    z.beta = huge

    return z
end
