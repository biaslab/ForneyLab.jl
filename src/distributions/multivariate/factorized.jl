export FactorizedDistribution

"""
Description:

    Encodes a distribution consisting of multiple factors.
    `p(x1,x2,x3) = p(x1)p(x2)p(x3)`
    The factors should have identical distribution types and dimensionality.

Parameters:

    factors::Vector{ProbabilityDistribution}

Construction:

    FactorizedDistribution([GaussianDistribution(), GaussianDistribution()])
"""
type FactorizedDistribution{dtype<:ProbabilityDistribution,n_factors} <: MultivariateProbabilityDistribution
    factors::Vector{dtype}

    function FactorizedDistribution(factors::Vector{dtype})
        (length(factors) > 1) || error("FactorizedDistribution should contain at least 2 factors")
        for i=2:length(factors)
            (typeof(factors[i]) == dtype) || error("All factors in FactorizedDistribution should have the same distribution type")
        end

        return new{dtype,length(factors)}(factors)
    end
end
FactorizedDistribution{T<:ProbabilityDistribution}(factors::Vector{T}) = FactorizedDistribution{typeof(factors[1]), length(factors)}(factors)
FactorizedDistribution() = FactorizedDistribution([GaussianDistribution(), GaussianDistribution()])

function show(io::IO, dist::FactorizedDistribution)
    println(io, typeof(dist))
    println(io, dist.factors)
end

function vague!{dtype,n_factors}(dist::FactorizedDistribution{dtype,n_factors})
    map(vague!, dist.factors)

    return dist
end

vague{dtype,n_factors}(::Type{FactorizedDistribution{dtype,n_factors}}) = FactorizedDistribution(dtype[vague(dtype) for i=1:n_factors])

Base.mean(dist::FactorizedDistribution) = vcat(map(mean, dist.factors)...)

sample(dist::FactorizedDistribution) = vcat(map(sample, dist.factors)...)

isProper(dist::FactorizedDistribution) = all(map(isProper, dist.factors))

function ==(x::FactorizedDistribution, y::FactorizedDistribution)
    (length(x.factors)==length(x.factors)) || return false
    return all(i->(x.factors[i]==y.factors[i]), collect(1:length(x.factors)))
end
