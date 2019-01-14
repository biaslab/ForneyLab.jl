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

ForneyLab.slug(::Type{Poisson}) = "Poisson"

ForneyLab.format(dist::ProbabilityDistribution{Univariate, Poisson}) = "$(slug(Poisson))(l=$(format(dist.params[:l])))"

ForneyLab.ProbabilityDistribution(::Type{Univariate}, ::Type{Poisson}; l=1.0) = ProbabilityDistribution{Univariate, Poisson}(Dict(:l=>l))
ForneyLab.ProbabilityDistribution(::Type{Poisson}; l=1.0) = ProbabilityDistribution{Univariate, Poisson}(Dict(:l=>l))

ForneyLab.dims(dist::ProbabilityDistribution{Univariate, Poisson}) = 1

ForneyLab.vague(::Type{Poisson}) = ProbabilityDistribution(Univariate, Poisson, l=huge)

ForneyLab.isProper(dist::ProbabilityDistribution{Univariate, Poisson}) = (0 < dist.params[:l] < huge)

ForneyLab.unsafeMean(dist::ProbabilityDistribution{Univariate, Poisson}) = dist.params[:l]

ForneyLab.unsafeVar(dist::ProbabilityDistribution{Univariate, Poisson}) = dist.params[:l]

ForneyLab.unsafeLogMean(dist::ProbabilityDistribution{Univariate, Poisson}) = dist.params[:l]

# TODO: skip the implementation of entropy and average energy for now, they require additional work

# Entropy functional
# function ForneyLab.differentialEntropy(dist::ProbabilityDistribution{Univariate, Poisson})
#     ...
# end

# Average energy functional
function ForneyLab.averageEnergy(::Type{Poisson}, marg_out::ProbabilityDistribution{Univariate}, marg_l::ProbabilityDistribution{Univariate})
    # approximation of expectation of logX!
    c_(j) = sum([binomial(j, k)*(-1)^(-k)*log(factorial(k)) for k in collect(0:j)])
    Elog_fact(λ, lim=20) = sum([(-λ)^j/factorial(j)*c_(j) for j in collect(0:lim)])

    -unsafeMean(marg_out)*unsafeLogMean(marg_l)+unsafeMean(marg_l)+Elog_fact(unsafeMean(marg_out))
end
