############################################
# InverseGammaDistribution
############################################
# Description:
#   Encodes an inverse gamma PDF.
#   Pamameters: scalars a (shape) and b (rate).
############################################
export InverseGammaDistribution

type InverseGammaDistribution <: ProbabilityDistribution
    a::Float64 # shape
    b::Float64 # rate
    InverseGammaDistribution(; a=1.0, b=1.0) = new(a, b)
end

vague(::Type{InverseGammaDistribution}) = InverseGammaDistribution(a=-1.0+tiny(), b=tiny())
function Base.mean(dist::InverseGammaDistribution)
    if dist.a > 1.0
        return dist.b / (dist.a - 1)
    else
        return NaN
    end
end
function Base.var(dist::InverseGammaDistribution)
    if dist.a > 2.0
        return (dist.b^2) / ( ( (dist.a-1)^2 ) * (dist.a-2) )
    else
        return NaN
    end
end

function show(io::IO, dist::InverseGammaDistribution)
    println(io, typeof(dist))
    println(io, "a = $(dist.a) (shape)")
    println(io, "b = $(dist.b) (rate)")
end

==(x::InverseGammaDistribution, y::InverseGammaDistribution) = (x.a==y.a && x.b==y.b)