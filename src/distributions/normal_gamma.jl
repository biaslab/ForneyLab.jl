############################################
# NormalGammaDistribution
############################################
# Description:
#   Encodes a normal-gamma PDF.
#   Pamameters: m (location), beta (precision) scalars a (shape) and b (rate).
############################################
export NormalGammaDistribution

type NormalGammaDistribution <: ProbabilityDistribution
    # All univariate, so parameters are floats
    m::Float64 # location
    beta::Float64 # precision
    a::Float64 # shape
    b::Float64 # rate
    NormalGammaDistribution(; m=0.0, beta=1.0, a=1.0, b=1.0) = new(m, beta, a, b)
end

vague(::Type{NormalGammaDistribution}) = NormalGammaDistribution(m=0.0, beta=1.0, a=1.0-tiny, b=tiny)

format(dist::NormalGammaDistribution) = "Ng(m=$(format(dist.m)), β=$(format(dist.beta)), a=$(format(dist.a)), b=$(format(dist.b)))"
show(io::IO, dist::NormalGammaDistribution) = println(io, format(dist))

==(x::NormalGammaDistribution, y::NormalGammaDistribution) = (x.m==y.m && x.beta==y.beta && x.a==y.a && x.b==y.b)