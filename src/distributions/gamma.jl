############################################
# GammaDistribution
############################################
# Description:
#   Encodes a gamma PDF.
#   Pamameters: scalars a (shape) and b (rate).
############################################
export GammaDistribution

type GammaDistribution <: ProbabilityDistribution
    a::Float64 # shape
    b::Float64 # rate
    GammaDistribution(; a=1.0, b=1.0) = new(a, b)
end

vague(::Type{GammaDistribution}) = GammaDistribution(a=tiny, b=tiny) # Scale invariant (Jeffrey's) prior

function Base.mean(dist::GammaDistribution)
    if dist.a > 0 && dist.b > 0
        return dist.a / dist.b
    else 
        return NaN
    end
end

function Base.var(dist::GammaDistribution)
    if dist.a > 0 && dist.b > 0
        return dist.a / (dist.b^2)
    else 
        return NaN
    end
end

format(dist::GammaDistribution) = "Gam(a=$(format(dist.a)), b=$(format(dist.b)))"
show(io::IO, dist::GammaDistribution) = println(io, format(dist))

==(x::GammaDistribution, y::GammaDistribution) = (x.a==y.a && x.b==y.b)