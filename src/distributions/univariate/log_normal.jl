export LogNormal

"""
Description:

    Encodes a log-normal PDF.

Pamameters:

    Real scalars m (location) and s=σ^2 (σ scale).

Construction:

    LogNormal(m=0.0, s=1.0)
"""
type LogNormal <: Univariate
    m::Float64 # location
    s::Float64 # squared scale (s=σ^2)
end

LogNormal(; m=0.0, s=1.0) = LogNormal(m, s)

function vague!(dist::LogNormal)
    dist.m = 0.0
    dist.s = huge
    return dist
end

isProper(dist::LogNormal) = dist.s >= tiny

Base.mean(dist::LogNormal) = isProper(dist) ? exp(dist.m + 0.5*dist.s) : NaN

Base.var(dist::LogNormal) = isProper(dist) ? (exp(dist.s) - 1)*exp(2*dist.m + dist.s) : NaN

format(dist::LogNormal) = "logN(μ=$(format(dist.m)), σ²=$(format(dist.s)))"

show(io::IO, dist::LogNormal) = println(io, format(dist))

==(x::LogNormal, y::LogNormal) = (x.m==y.m && x.s==y.s)

@symmetrical function prod!(x::LogNormal, y::Delta{Float64}, z::Delta{Float64}=Delta(y.m))
    # Product of log-normal PDF and Delta
    (y.m >= 0.) || throw(DomainError())
    (z.m == y.m) || (z.m = y.m)

    return z
end

@symmetrical function prod!(x::LogNormal, y::Delta{Float64}, z::LogNormal)
    # Product of multivariate log-normal PDF and MvDelta, force result to be mv log-normal
    (y.m >= 0.) || throw(DomainError())
    z.m = log(y.m)
    z.s = tiny

    return z
end
