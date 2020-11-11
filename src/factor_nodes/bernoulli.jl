export Bernoulli

"""
Description:

    Bernoulli factor node

    out ∈ {0, 1}
    p ∈ [0, 1]

    f(out, p) = Ber(out|p) = p^out (1 - p)^{1 - out}

Interfaces:

    1. out
    2. p

Construction:

    Bernoulli(id=:some_id)
"""
mutable struct Bernoulli <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function Bernoulli(out, p; id=generateId(Bernoulli))
        @ensureVariables(out, p)
        self = new(id, Array{Interface}(undef, 2), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:p] = self.interfaces[2] = associate!(Interface(self), p)

        return self
    end
end

slug(::Type{Bernoulli}) = "Ber"

format(dist::ProbabilityDistribution{Univariate, Bernoulli}) = "$(slug(Bernoulli))(p=$(format(dist.params[:p])))"

ProbabilityDistribution(::Type{Univariate}, ::Type{Bernoulli}; p=0.5) = ProbabilityDistribution{Univariate, Bernoulli}(Dict(:p=>p))
ProbabilityDistribution(::Type{Bernoulli}; p=0.5) = ProbabilityDistribution{Univariate, Bernoulli}(Dict(:p=>p))

dims(dist::ProbabilityDistribution{Univariate, Bernoulli}) = 1

vague(::Type{Bernoulli}) = ProbabilityDistribution(Univariate, Bernoulli, p=0.5)

isProper(dist::ProbabilityDistribution{Univariate, Bernoulli}) = (0 <= dist.params[:p] <= 1)

unsafeMean(dist::ProbabilityDistribution{Univariate, Bernoulli}) = dist.params[:p]

unsafeMode(dist::ProbabilityDistribution{Univariate, Bernoulli}) = round(dist.params[:p])

unsafeMeanVector(dist::ProbabilityDistribution{Univariate, Bernoulli}) = [dist.params[:p], 1 - dist.params[:p]]

unsafeVar(dist::ProbabilityDistribution{Univariate, Bernoulli}) = dist.params[:p]*(1 - dist.params[:p])

logPdf(dist::ProbabilityDistribution{Univariate, Bernoulli}, x) = x*log(dist.params[:p]) + (1.0-x)*log(1.0-dist.params[:p])

sample(dist::ProbabilityDistribution{Univariate, Bernoulli}) = 1.0*(rand() < dist.params[:p])

function prod!( x::ProbabilityDistribution{Univariate, Bernoulli},
                y::ProbabilityDistribution{Univariate, Bernoulli},
                z::ProbabilityDistribution{Univariate, Bernoulli}=ProbabilityDistribution(Univariate, Bernoulli, p=0.5))

    norm = x.params[:p] * y.params[:p] + (1 - x.params[:p]) * (1 - y.params[:p])
    (norm > 0) || error("Product of $(x) and $(y) cannot be normalized")
    z.params[:p] = (x.params[:p] * y.params[:p]) / norm

    return z
end

# Entropy functional
function differentialEntropy(dist::ProbabilityDistribution{Univariate, Bernoulli})
    p = clamp(dist.params[:p], tiny, 1.0 - tiny)

    -(1.0 - p)*log(1.0 - p) -
    p*log(p)
end

# Average energy functional
function averageEnergy(::Type{Bernoulli}, marg_out::ProbabilityDistribution{Univariate}, marg_p::ProbabilityDistribution{Univariate})
    -unsafeMean(marg_out)*unsafeLogMean(marg_p) -
    (1.0 - unsafeMean(marg_out))*unsafeMirroredLogMean(marg_p)
end
