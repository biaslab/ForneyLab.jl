export Bernoulli

"""
Description:

    Encodes a distribution over binary domain {false,true}.

Pamameters:

    p ∈ [0,1]: Pr{X=true} = p

Construction:

    Bernoulli(p)
"""
type Bernoulli <: Univariate
    p::Float64 # Pr{X=true}
end

Bernoulli() = Bernoulli(0.5)

pdf(dist::Bernoulli, x::Bool) = (x ? dist.p : (1.0-dist.p))

function vague!(dist::Bernoulli)
    dist.p = 0.5
    return dist
end

isProper(dist::Bernoulli) = (0 <= dist.p <= 1)

unsafeMean(dist::Bernoulli) = dist.p

unsafeVar(dist::Bernoulli) = dist.p*(1-dist.p)

sample(dist::Bernoulli) = (rand() < dist.p)

format(dist::Bernoulli) = "Bernoulli(p=$(format(dist.p)))"

show(io::IO, dist::Bernoulli) = println(io, format(dist))

==(x::Bernoulli, y::Bernoulli) = isApproxEqual(x.p, y.p)

function prod!(x::Bernoulli, y::Bernoulli, z::Bernoulli=Bernoulli())
    # Multiplication of 2 Bernoulli PDFs: p(z) = p(x) * p(y)
    norm = x.p * y.p + (1 - x.p) * (1 - y.p)
    (norm > 0) || error("Product of $(x) and $(y) cannot be normalized")
    z.p = (x.p * y.p) / norm

    return z
end

@symmetrical function prod!(x::Bernoulli, y::Delta{Bool}, z::Delta{Bool}=Delta(true))
    # Product of Bernoulli PMF and Delta
    if y.m
        (x.p > 0.) || error("Invalid product of Bernoulli and Delta")
    else
        (x.p < 1.) || error("Invalid product of Bernoulli and Delta")
    end
    z.m = y.m

    return z
end

@symmetrical function prod!(x::Bernoulli, y::Delta{Bool}, z::Bernoulli)
    # Product of Bernoulli PMF and Delta, force result to be Bernoulli
    if y.m
        (x.p > 0.) || error("Invalid product of Bernoulli and Delta")
    else
        (x.p < 1.) || error("Invalid product of Bernoulli and Delta")
    end
    z.p = (y.m) ? 1.0 : 0.0

    return z
end

# Converts from Delta{Bool} -> Bernoulli
Base.convert{T<:Bool}(::Type{Bernoulli}, delta::Delta{T}) = Bernoulli(float(delta.m))

function differentialEntropy(dist::Bernoulli)
    -(1.0 - dist.p)*log(1.0 - dist.p) -
    dist.p*log(dist.p)
end