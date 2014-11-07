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

uninformative(::Type{NormalGammaDistribution}) = NormalGammaDistribution(m=0.0, beta=1.0, a=1.0-tiny(), b=tiny())

function show(io::IO, dist::NormalGammaDistribution)
    println(io, typeof(dist))
    println(io, "m = $(dist.m) (location)")
    println(io, "beta = $(dist.beta) (precision)")
    println(io, "a = $(dist.a) (shape)")
    println(io, "b = $(dist.b) (rate)")
end

==(x::NormalGammaDistribution, y::NormalGammaDistribution) = (x.m==y.m && x.beta==y.beta && x.a==y.a && x.b==y.b)