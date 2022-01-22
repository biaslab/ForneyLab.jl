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

    function Poisson(out, l; id=generateId(Poisson))
        @ensureVariables(out, l)
        self = new(id, Array{Interface}(undef, 2), Dict{Symbol, Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:l] = self.interfaces[2] = associate!(Interface(self), l)
        return self
    end
end

slug(::Type{Poisson}) = "Poisson"

format(dist::Distribution{Univariate, Poisson}) = "$(slug(Poisson))(l=$(format(dist.params[:l])))"

Distribution(::Type{Univariate}, ::Type{Poisson}; l=1.0) = Distribution{Univariate, Poisson}(Dict(:l=>l))
Distribution(::Type{Poisson}; l=1.0) = Distribution{Univariate, Poisson}(Dict(:l=>l))

dims(dist::Distribution{Univariate, Poisson}) = ()

vague(::Type{Poisson}) = Distribution(Univariate, Poisson, l=huge)

isProper(dist::Distribution{Univariate, Poisson}) = (0 < dist.params[:l] < huge)

unsafeMean(dist::Distribution{Univariate, Poisson}) = Float64(dist.params[:l])

unsafeVar(dist::Distribution{Univariate, Poisson}) = Float64(dist.params[:l])

logPdf(dist::Distribution{Univariate, Poisson}, x) = x*log(dist.params[:l]) - dist.params[:l] - logfactorial(x)

sample(dist::Distribution{Univariate, Poisson}) = poisinvcdf(dist.params[:l], rand())

naturalParams(dist::Distribution{Univariate, Poisson}) = [log(dist.params[:l])]

standardDistribution(V::Type{Univariate}, F::Type{Poisson}; η::Vector) = Distribution(V, F, l=exp(η[1]))

logNormalizer(::Type{Univariate}, ::Type{Poisson}; η::Vector) = exp(η[1])

logPdf(V::Type{Univariate}, F::Type{Poisson}, x::Number; η::Vector) = -logfactorial(x) + [x]'*η - logNormalizer(V, F; η=η)

# ∑ [λ^k*log(k!)]/k! from k=0 to inf
# Approximates the above sum for calculation of averageEnergy and differentialEntropy
# @ref https://arxiv.org/pdf/1708.06394.pdf
function apprSum(l, j=100)
    sum([(l)^(k)*logfactorial(k)/exp(logfactorial(k)) for k in collect(0:j)])
end

# Entropy functional
# @ref https://en.wikipedia.org/wiki/Poisson_distribution
function differentialEntropy(dist::Distribution{Univariate, Poisson})
    l = clamp(dist.params[:l], tiny, huge)
    l*(1-log(l)) + exp(-l)*apprSum(l)
end

# Average energy functional
function averageEnergy(::Type{Poisson}, marg_out::Distribution{Univariate}, marg_l::Distribution{Univariate})
    unsafeMean(marg_l) -
    unsafeMean(marg_out)*unsafeLogMean(marg_l) +
    exp(-unsafeMean(marg_out))*apprSum(unsafeMean(marg_out))
end
