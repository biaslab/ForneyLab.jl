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

isProper(dist::LogNormalDistribution) = dist.s > tiny

Base.mean(dist::LogNormalDistribution) = isProper(dist) ? exp(dist.m + 0.5*dist.s) : NaN

Base.mean(::Type{DeltaDistribution{Float64}}, d::LogNormalDistribution) = DeltaDistribution(mean(d)) # Definition for post-processing

Base.var(dist::LogNormalDistribution) = isProper(dist) ? (exp(dist.s) - 1)*exp(2*dist.m + dist.s) : NaN

format(dist::LogNormalDistribution) = "logN(μ=$(format(dist.m)), σ²=$(format(dist.s)))"

show(io::IO, dist::LogNormalDistribution) = println(io, format(dist))

==(x::LogNormalDistribution, y::LogNormalDistribution) = (x.m==y.m && x.s==y.s)
