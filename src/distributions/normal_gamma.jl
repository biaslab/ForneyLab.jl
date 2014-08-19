############################################
# NormalGammaDistribution
############################################
# Description:
#   Encodes a normal-gamma PDF.
#   Pamameters: m (location), W (precision) scalars a (shape) and b (rate).
############################################
export NormalGammaDistribution

type NormalGammaDistribution <: ProbabilityDistribution
    # All univariate, so parameters are floats
    m::Float64 # location
    W::Float64 # precision
    a::Float64 # shape
    b::Float64 # rate
    NormalGammaDistribution(; m=0.0, W=1.0, a=1.0, b=1.0) = new(m, W, a, b)
end

uninformative(dist_type::Type{NormalGammaDistribution}) = NormalGammaDistribution(m=0.0, W=1.0, a=0.999, b=0.001)

function show(io::IO, dist::NormalGammaDistribution)
    println(io, typeof(dist))
    println(io, "m = $(dist.m) (location)")
    println(io, "W = $(dist.W) (precision)")
    println(io, "a = $(dist.a) (shape)")
    println(io, "b = $(dist.b) (rate)")
end

==(x::NormalGammaDistribution, y::NormalGammaDistribution) = (x.m==y.m && x.W==y.W && x.a==y.a && x.b==y.b)