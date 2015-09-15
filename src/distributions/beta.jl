############################################
# BetaDistribution
############################################
# Description:
#   Encodes a beta PDF.
#   Pamameters: a and b.
############################################
export BetaDistribution

type BetaDistribution <: ProbabilityDistribution
    a::Float64 # shape
    b::Float64 # rate
    BetaDistribution(; a=1.0, b=1.0) = new(a, b)
end

vague(::Type{BetaDistribution}) = BetaDistribution(a=tiny, b=tiny)
function Base.mean(dist::BetaDistribution)
    (a > 0 && b > 0) || return NaN
    return a/(a+b)
end
function Base.var(dist::BetaDistribution)
    (a > 0 && b > 0) || return NaN
    return a*b/((a+b)^2*(a+b+1))
end

format(dist::BetaDistribution) = "Bet(a=$(format(dist.a)), b=$(format(dist.b)))"
show(io::IO, dist::BetaDistribution) = println(io, format(dist))

==(x::BetaDistribution, y::BetaDistribution) = (x.a==y.a && x.b==y.b)