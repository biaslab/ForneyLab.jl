export BernoulliDistribution

"""
Description:

    Encodes a distribution over binary domain {false,true}.

Pamameters:

    p ∈ [0,1]: Pr{X=true} = p

Construction:

    BernoulliDistribution(p)
"""
type BernoulliDistribution <: UnivariateProbabilityDistribution
    p::Float64 # Pr{X=true}
end
BernoulliDistribution() = BernoulliDistribution(0.5)

function vague!(dist::BernoulliDistribution)
    dist.p = 0.5
    return dist
end

isProper(dist::BernoulliDistribution) = (0 <= dist.p <= 1)

Base.mean(dist::BernoulliDistribution) = dist.p

Base.var(dist::BernoulliDistribution) = dist.p*(1-dist.p)

sample(dist::BernoulliDistribution) = (rand() < dist.p)

format(dist::BernoulliDistribution) = "Bernoulli(p=$(format(dist.p)))"

show(io::IO, dist::BernoulliDistribution) = println(io, format(dist))

==(x::BernoulliDistribution, y::BernoulliDistribution) = isApproxEqual(x.p, y.p)

function prod!(x::BernoulliDistribution, y::BernoulliDistribution, z::BernoulliDistribution=BernoulliDistribution())
    # Multiplication of 2 Bernoulli PDFs: p(z) = p(x) * p(y)
    norm = x.p * y.p + (1 - x.p) * (1 - y.p)
    (norm > 0) || error("Product of $(x) and $(y) cannot be normalized")
    z.p = (x.p * y.p) / norm

    return z
end

@symmetrical function prod!(x::BernoulliDistribution, y::DeltaDistribution{Bool}, z::DeltaDistribution{Bool}=DeltaDistribution(true))
    # Product of Bernoulli PMF and Delta
    if y.m
        (x.p > 0.) || error("Invalid product of Bernoulli and Delta")
    else
        (x.p < 1.) || error("Invalid product of Bernoulli and Delta")
    end
    z.m = y.m

    return z
end

@symmetrical function prod!(x::BernoulliDistribution, y::DeltaDistribution{Bool}, z::BernoulliDistribution)
    # Product of Bernoulli PMF and Delta, force result to be Bernoulli
    if y.m
        (x.p > 0.) || error("Invalid product of Bernoulli and Delta")
    else
        (x.p < 1.) || error("Invalid product of Bernoulli and Delta")
    end
    z.p = (y.m) ? 1.0 : 0.0

    return z
end

# Converts from DeltaDistribution{Bool} -> BernoulliDistribution
Base.convert{T<:Bool}(::Type{BernoulliDistribution}, delta::DeltaDistribution{T}) = BernoulliDistribution(float(delta.m))
