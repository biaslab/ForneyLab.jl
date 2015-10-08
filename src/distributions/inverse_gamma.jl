############################################
# InverseGammaDistribution
############################################
# Description:
#   Encodes an inverse gamma PDF.
#   Pamameters: scalars a (shape) and b (scale).
############################################
export InverseGammaDistribution

type InverseGammaDistribution <: ProbabilityDistribution
    a::Float64 # shape
    b::Float64 # scale
    InverseGammaDistribution(; a=3.0, b=2.0) = new(a, b)
end

vague(::Type{InverseGammaDistribution}) = InverseGammaDistribution(a=tiny, b=tiny) # Jeffrey's prior
isProper(dist::InverseGammaDistribution) = (dist.a >= tiny && dist.b >= tiny)
function Base.mean(dist::InverseGammaDistribution)
    if isProper(dist) && dist.a > 1.0
        return dist.b / (dist.a - 1)
    else
        return NaN
    end
end
function Base.var(dist::InverseGammaDistribution)
    if isProper(dist) && dist.a > 2.0
        return (dist.b^2) / ( ( (dist.a-1)^2 ) * (dist.a-2) )
    else
        return NaN
    end
end

format(dist::InverseGammaDistribution) = "Ig(a=$(format(dist.a)), b=$(format(dist.b)))"
show(io::IO, dist::InverseGammaDistribution) = println(io, format(dist))

==(x::InverseGammaDistribution, y::InverseGammaDistribution) = (x.a==y.a && x.b==y.b)