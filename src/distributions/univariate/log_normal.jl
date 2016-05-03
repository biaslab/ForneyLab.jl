export LogNormalDistribution

"""
Description:

    Encodes a log-normal PDF.

Pamameters:

    Real scalars m (location) and s=σ^2 (σ scale).

Construction:

    LogNormalDistribution(m=0.0, s=1.0)
"""
type LogNormalDistribution <: UnivariateProbabilityDistribution
    m::Float64 # location
    s::Float64 # squared scale (s=σ^2)
end

LogNormalDistribution(; m=0.0, s=1.0) = LogNormalDistribution(m, s)

function vague!(dist::LogNormalDistribution)
    dist.m = 0.0
    dist.s = huge
    return dist
end

isProper(dist::LogNormalDistribution) = dist.s >= tiny

Base.mean(dist::LogNormalDistribution) = isProper(dist) ? exp(dist.m + 0.5*dist.s) : NaN

Base.var(dist::LogNormalDistribution) = isProper(dist) ? (exp(dist.s) - 1)*exp(2*dist.m + dist.s) : NaN

format(dist::LogNormalDistribution) = "logN(μ=$(format(dist.m)), σ²=$(format(dist.s)))"

show(io::IO, dist::LogNormalDistribution) = println(io, format(dist))

==(x::LogNormalDistribution, y::LogNormalDistribution) = (x.m==y.m && x.s==y.s)

@symmetrical function prod!(x::LogNormalDistribution, y::DeltaDistribution{Float64}, z::DeltaDistribution{Float64}=DeltaDistribution(y.m))
    # Product of log-normal PDF and Delta
    (y.m >= 0.) || throw(DomainError())
    (z.m == y.m) || (z.m = y.m)

    return z
end

@symmetrical function prod!(x::LogNormalDistribution, y::DeltaDistribution{Float64}, z::LogNormalDistribution)
    # Product of multivariate log-normal PDF and MvDelta, force result to be mv log-normal
    (y.m >= 0.) || throw(DomainError())
    z.m = log(y.m)
    z.s = tiny

    return z
end
