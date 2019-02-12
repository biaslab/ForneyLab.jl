export Poisson
import SpecialFunctions.lfactorial

"""
Description:
    Poisson factor node

    Real scalars
    l > 0 (rate)

    f(out, l) = Poisson(out|l) = 1/(x!) * l^x * exp(-l)

Interfaces:
    1. out
    2. l

Construction:
    Poisson(id=:some_id)
"""
mutable struct Poisson <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function Poisson(out, l, id=generateId(Poisson))
        @ensureVariables(out, l)
        self = new(id, Array{Interface}(undef, 2), Dict{Symbol, Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:l] = self.interfaces[2] = associate!(Interface(self), l)
        return self
    end
end

slug(::Type{Poisson}) = "Poisson"

format(dist::ProbabilityDistribution{Univariate, Poisson}) = "$(slug(Poisson))(l=$(format(dist.params[:l])))"

ProbabilityDistribution(::Type{Univariate}, ::Type{Poisson}; l=1.0) = ProbabilityDistribution{Univariate, Poisson}(Dict(:l=>l))
ProbabilityDistribution(::Type{Poisson}; l=1.0) = ProbabilityDistribution{Univariate, Poisson}(Dict(:l=>l))

dims(dist::ProbabilityDistribution{Univariate, Poisson}) = 1

vague(::Type{Poisson}) = ProbabilityDistribution(Univariate, Poisson, l=huge)

isProper(dist::ProbabilityDistribution{Univariate, Poisson}) = (0 < dist.params[:l] < huge)

unsafeMean(dist::ProbabilityDistribution{Univariate, Poisson}) = dist.params[:l]

unsafeVar(dist::ProbabilityDistribution{Univariate, Poisson}) = dist.params[:l]

# ∑ [λ^k*log(k!)]/k! from k=0 to inf
# Approximates the above sum for differentialEntropy calculation
apprSum(l, j=66) = sum([(l)^(k)*lfactorial(k)/exp(lfactorial(k)) for k in collect(0:j)])

# Entropy functional
# @ref https://en.wikipedia.org/wiki/Poisson_distribution
function differentialEntropy(dist::ProbabilityDistribution{Univariate, Poisson})

    l = clamp(dist.params[:l], tiny, huge)
    l*(1-log(l)) + exp(-l)*apprSum(l)
end


# ∑ binomial(j, k)*(-1)^{-k}*log(k!) k=0 to inf
# logarithmic difference coefficient
logDiffCoef(j::Int64) = sum([binomial(j, k)*(-1)^(-k)*lfactorial(k) for k in collect(0:j)])
# approximation of expectation of logX!
# @ref https://arxiv.org/pdf/1708.06394.pdf
unsafeLogFact(l, lim=66) = sum([(-l)^j/exp(lfactorial(j))*logDiffCoef(j) for j in collect(0:lim)])

# Average energy functional
function averageEnergy(::Type{Poisson}, marg_out::ProbabilityDistribution{Univariate}, marg_l::ProbabilityDistribution{Univariate})

    -unsafeMean(marg_out)*unsafeLogMean(marg_l)+unsafeMean(marg_l)+unsafeLogFact(unsafeMean(marg_out))
end
