export InverseGamma

"""
Description:

    Encodes an inverse gamma PDF.

Pamameters:

    Real scalars a > 0 (shape) and b > 0 (scale).

Construction:

    InverseGamma(a=1.0, b=1.0)

Reference:

    Korl, 2005; A factor graph approach to signal modelling, system identification and filtering; appendix A
"""
type InverseGamma <: Univariate
    a::Float64 # shape
    b::Float64 # scale
end

InverseGamma(; a=3.0, b=2.0) = InverseGamma(a, b)

function pdf(dist::InverseGamma, x::Float64)
    (0.0 <= x) || return 0.0
    C = (dist.b^dist.a) / gamma(dist.a)
    return C * x^(-dist.a-1) * exp(-dist.b/x)
end

function vague!(dist::InverseGamma)
    dist.a = tiny
    dist.b = tiny
    return dist
end

isProper(dist::InverseGamma) = (dist.a >= tiny && dist.b >= tiny)

function Base.mean(dist::InverseGamma)
    if isProper(dist) && dist.a > 1.0
        return dist.b / (dist.a - 1)
    else
        return NaN
    end
end

function Base.var(dist::InverseGamma)
    if isProper(dist) && dist.a > 2.0
        return (dist.b^2) / ( ( (dist.a-1)^2 ) * (dist.a-2) )
    else
        return NaN
    end
end

format(dist::InverseGamma) = "Ig(a=$(format(dist.a)), b=$(format(dist.b)))"

show(io::IO, dist::InverseGamma) = println(io, format(dist))

function prod!(x::InverseGamma, y::InverseGamma, z::InverseGamma=InverseGamma())
    # Multiplication of 2 inverse gamma PDFs: p(z) = p(x) * p(y)
    z.a = x.a + y.a + 1.0
    z.b = x.b + y.b

    return z
end

@symmetrical function prod!(x::InverseGamma, y::Delta{Float64}, z::Delta{Float64}=Delta(1.0))
    # Multiplication of inverse gamma PDF with a Delta
    (y.m >= 0.0) || throw(DomainError())
    z.m = y.m

    return z
end

@symmetrical function prod!(x::InverseGamma, y::Delta{Float64}, z::InverseGamma)
    # Multiplication of inverse gamma PDF with a Delta, force result to be inverse gamma
    (y.m >= 0.0) || throw(DomainError())
    z.a = clamp(1e8*y.m, tiny, huge)
    z.b = clamp(y.m*(z.a-1), tiny, huge)

    return z
end

==(x::InverseGamma, y::InverseGamma) = (x.a==y.a && x.b==y.b)