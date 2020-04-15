export Poisson

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

unsafeMean(dist::ProbabilityDistribution{Univariate, Poisson}) = Float64(dist.params[:l])

unsafeVar(dist::ProbabilityDistribution{Univariate, Poisson}) = Float64(dist.params[:l])

logPdf(dist::ProbabilityDistribution{Univariate, Poisson}, x) = x*log(dist.params[:l]) - dist.params[:l] - logfactorial(x)

sample(dist::ProbabilityDistribution{Univariate, Poisson}) = poisinvcdf(dist.params[:l], rand())

# ∑ [λ^k*log(k!)]/k! from k=0 to inf
# Approximates the above sum for calculation of averageEnergy and differentialEntropy
# @ref https://arxiv.org/pdf/1708.06394.pdf
function apprSum(l, j=100)
    sum([(l)^(k)*logfactorial(k)/exp(logfactorial(k)) for k in collect(0:j)])
end

# Entropy functional
# @ref https://en.wikipedia.org/wiki/Poisson_distribution
function differentialEntropy(dist::ProbabilityDistribution{Univariate, Poisson})
    l = clamp(dist.params[:l], tiny, huge)
    l*(1-log(l)) + exp(-l)*apprSum(l)
end

# Average energy functional
function averageEnergy(::Type{Poisson}, marg_out::ProbabilityDistribution{Univariate}, marg_l::ProbabilityDistribution{Univariate})
    unsafeMean(marg_l) -
    unsafeMean(marg_out)*unsafeLogMean(marg_l) +
    exp(-unsafeMean(marg_out))*apprSum(unsafeMean(marg_out))
end
