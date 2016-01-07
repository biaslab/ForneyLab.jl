############################################
# LogNormalDistribution
############################################
# Description:
#   Encodes a log-normal PDF.
#   Pamameters: scalars m (location) and s=σ^2 (σ scale).
############################################
export LogNormalDistribution, approximateWithGamma, approximateWithLogNormal

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

# Post-processing functions for approximations
# Approximations are computed with moment matching
function approximateWithGamma(d::LogNormalDistribution)
    a = 1/(exp(d.s) - 1)
    b = exp(-(d.m + 0.5*d.s))/(exp(d.s)-1)
    return GammaDistribution(a=a, b=b)
end

function approximateWithLogNormal(d::GammaDistribution)
    m = (3/2)*log(d.a) - 0.5*log(1 + d.a) - log(d.b)
    s = log(1 + d.a) - log(d.a)
    return LogNormalDistribution(m=m, s=s)
end