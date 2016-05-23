export Beta

"""
Description:

    Encodes a beta PDF.

Pamameters:

    Real scalars a > 0 (shape) and b > 0 (rate).

Construction:

    Beta(a=1.0, b=1.0)

Reference:

    Bishop, 2006; Pattern recognition and machine learning; appendix B
"""
type Beta <: Univariate
    a::Float64 # shape
    b::Float64 # rate
end

Beta(; a=1.0, b=1.0) = Beta(a, b)

function pdf(dist::Beta, x::Float64)
    (0.0 <= x <= 1.0) || return 0.0
    B = gamma(dist.a)*gamma(dist.b) / gamma(dist.a+dist.b)
    return (x^(dist.a-1)*(1-x)^(dist.b-1)) / B
end


function vague!(dist::Beta)
    dist.a = tiny
    dist.b = tiny
    return dist
end

isProper(dist::Beta) = (dist.a >= tiny && dist.b >= tiny)

Base.mean(dist::Beta) = isProper(dist) ? dist.a/(dist.a+dist.b) : NaN

Base.var(dist::Beta) = isProper(dist) ? dist.a*dist.b/((dist.a+dist.b)^2*(dist.a+dist.b+1)) : NaN

format(dist::Beta) = "Bet(a=$(format(dist.a)), b=$(format(dist.b)))"

show(io::IO, dist::Beta) = println(io, format(dist))

function prod!(x::Beta, y::Beta, z::Beta=Beta())
    # Multiplication of 2 beta PDFs: p(z) = p(x) * p(y)
    z.a = x.a + y.a - 1.0
    z.b = x.b + y.b - 1.0

    return z
end

@symmetrical function prod!(x::Beta, y::Delta{Float64}, z::Delta{Float64}=Delta(1.0))
    # Multiplication of beta PDF with Delta
    (0.0 <= y.m <= 1.0) || throw(DomainError())
    z.m = y.m

    return z
end

@symmetrical function prod!(x::Beta, y::Delta{Float64}, z::Beta)
    # Multiplication of beta PDF with Delta, force result to be beta
    (0.0 <= y.m <= 1.0) || throw(DomainError())
    z.a = (1.0 - y.m) * huge
    z.b = y.m * huge

    return z
end

==(x::Beta, y::Beta) = (x.a==y.a && x.b==y.b)
