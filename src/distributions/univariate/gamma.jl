export Gamma

"""
Description:

    Encodes a gamma PDF.

Pamameters:

    Real scalars scalars a > 0 (shape) and b > 0 (rate).

Construction:

    Gamma(a=1.0, b=1.0)

Reference:

    Bishop, 2006; Pattern recognition and machine learning; appendix B
"""
type Gamma <: Univariate
    a::Float64 # shape
    b::Float64 # rate
end

Gamma(; a=1.0, b=1.0) = Gamma(a, b)

function pdf(dist::Gamma, x::Float64)
    (0.0 <= x) || return 0.0
    C = (dist.b^dist.a) / gamma(dist.a)
    return C * x^(dist.a-1) * exp(-dist.b*x)
end

function vague!(dist::Gamma)
    dist.a = tiny
    dist.b = tiny
    return dist
end

isProper(dist::Gamma) = (dist.a >= tiny && dist.b >= tiny)

m(dist::Gamma) = dist.a/dist.b

V(dist::Gamma) = dist.a/(dist.b^2)

format(dist::Gamma) = "Gam(a=$(format(dist.a)), b=$(format(dist.b)))"

show(io::IO, dist::Gamma) = println(io, format(dist))

function prod!(x::Gamma, y::Gamma, z::Gamma=Gamma())
    # Multiplication of 2 gamma PDFs: p(z) = p(x) * p(y)
    z.a = x.a + y.a - 1.0
    z.b = x.b + y.b

    return z
end

@symmetrical function prod!(x::Gamma, y::Delta{Float64}, z::Delta{Float64}=Delta(1.0))
    # Multiplication of Gamma PDF with a Delta
    (y.m >= 0.0) || throw(DomainError())
    z.m = y.m

    return z
end

@symmetrical function prod!(x::Gamma, y::Delta{Float64}, z::Gamma)
    # Multiplication of Gamma PDF with a Delta, force result to be Gamma
    (y.m >= 0.0) || throw(DomainError())
    z.b = clamp(y.m / 1e-10, tiny, huge)
    z.a = clamp(y.m * z.b, tiny, huge)

    return z
end

@symmetrical function prod!(::Void, y::Delta{Float64}, z::Gamma)
    # Multiplication of an unknown with Delta, force result to be Gamma
    (y.m >= 0.0) || throw(DomainError())
    z.b = clamp(y.m / 1e-10, tiny, huge)
    z.a = clamp(y.m * z.b, tiny, huge)

    return z
end

==(x::Gamma, y::Gamma) = (x.a==y.a && x.b==y.b)

# Entropy functional
function H(dist::Gamma)
    return  log(gamma(dist.a)) -
            (dist.a - 1.0)*digamma(dist.a) -
            log(dist.b) +
            dist.a
end