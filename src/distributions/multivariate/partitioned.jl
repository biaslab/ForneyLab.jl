export PartitionedDistribution

"""
Description:

    Encodes the product of multiple independent distributions of the same type.
    `p(x1,x2,x3) = p(x1)p(x2)p(x3)`

Parameters:

    dtype<:ProbabilityDistribution: type of the factors
    n_factors: number of factors

Fields:

    factors::Vector{ProbabilityDistribution}

Construction:

    PartitionedDistribution([GaussianDistribution(), GaussianDistribution()])
"""
type PartitionedDistribution{dtype<:ProbabilityDistribution,n_factors} <: MultivariateProbabilityDistribution
    factors::Vector{dtype}

    function PartitionedDistribution(factors::Vector{dtype})
        (length(factors) > 1) || error("PartitionedDistribution should contain at least 2 factors")
        for i=2:length(factors)
            (typeof(factors[i]) == dtype) || error("All factors in PartitionedDistribution should have the same type")
        end

        return new{dtype,length(factors)}(factors)
    end
end
PartitionedDistribution{T<:ProbabilityDistribution}(factors::Vector{T}) = PartitionedDistribution{typeof(factors[1]), length(factors)}(factors)
PartitionedDistribution() = PartitionedDistribution([GaussianDistribution(), GaussianDistribution()])

function show(io::IO, dist::PartitionedDistribution)
    println(io, typeof(dist))
    println(io, dist.factors)
end

function vague!{dtype,n_factors}(dist::PartitionedDistribution{dtype,n_factors})
    map(vague!, dist.factors)

    return dist
end

vague{dtype,n_factors}(::Type{PartitionedDistribution{dtype,n_factors}}) = PartitionedDistribution(dtype[vague(dtype) for i=1:n_factors])

Base.mean(dist::PartitionedDistribution) = vcat(map(mean, dist.factors)...)

sample(dist::PartitionedDistribution) = vcat(map(sample, dist.factors)...)

isProper(dist::PartitionedDistribution) = all(map(isProper, dist.factors))

function ==(x::PartitionedDistribution, y::PartitionedDistribution)
    (length(x.factors)==length(x.factors)) || return false
    return all(i->(x.factors[i]==y.factors[i]), collect(1:length(x.factors)))
end
