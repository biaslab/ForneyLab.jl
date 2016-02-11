export NormalGammaDistribution

"""
Description:

    Encodes a normal-gamma distribution.
    Pamameters: m (location), beta (precision) scalars a (shape) and b (rate).
"""
type NormalGammaDistribution <: MultivariateProbabilityDistribution
    # All univariate, so parameters are floats
    m::Float64    # location
    beta::Float64 # precision
    a::Float64    # shape
    b::Float64    # rate
end

NormalGammaDistribution(; m=0.0, beta=1.0, a=1.0, b=1.0) = NormalGammaDistribution(m, beta, a, b)

# TODO: reference
function vague!(dist::NormalGammaDistribution)
    dist.m = 0.0
    dist.beta = 1.0
    dist.a = tiny
    dist.b = tiny
    return dist
end

vague(::Type{NormalGammaDistribution}) = NormalGammaDistribution(m=0.0, beta=1.0, a=tiny, b=tiny)

isProper(dist::NormalGammaDistribution) = ((dist.beta >= tiny) && (dist.a >= tiny)  && (dist.b >= tiny))

format(dist::NormalGammaDistribution) = "Ng(m=$(format(dist.m)), β=$(format(dist.beta)), a=$(format(dist.a)), b=$(format(dist.b)))"

show(io::IO, dist::NormalGammaDistribution) = println(io, format(dist))

==(x::NormalGammaDistribution, y::NormalGammaDistribution) = (x.m==y.m && x.beta==y.beta && x.a==y.a && x.b==y.b)

dimensions(message::Message{NormalGammaDistribution}) = 2
dimensions(distribution::NormalGammaDistribution) = 2
dimensions(message_type::Type{Message{NormalGammaDistribution}}) = 2
dimensions(distribution_type::Type{NormalGammaDistribution}) = 2
