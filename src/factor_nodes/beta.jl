export Beta

"""
Description:

    Beta factor node

    Real scalars
    a > 0
    b > 0

    f(out, a, b) = Beta(out|a, b) = Γ(a + b)/(Γ(a) Γ(b)) out^{a - 1} (1 - out)^{b - 1}

Interfaces:

    1. out
    2. a
    3. b

Construction:

    Beta(id=:some_id)
"""
mutable struct Beta <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function Beta(out, a, b; id=generateId(Beta))
        @ensureVariables(out, a, b)
        self = new(id, Array{Interface}(undef, 3), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:a] = self.interfaces[2] = associate!(Interface(self), a)
        self.i[:b] = self.interfaces[3] = associate!(Interface(self), b)

        return self
    end
end

slug(::Type{Beta}) = "Beta"

format(dist::Distribution{Univariate, Beta}) = "$(slug(Beta))(a=$(format(dist.params[:a])), b=$(format(dist.params[:b])))"

Distribution(::Type{Univariate}, ::Type{Beta}; a=1.0, b=1.0) = Distribution{Univariate, Beta}(Dict(:a=>a, :b=>b))
Distribution(::Type{Beta}; a=1.0, b=1.0) = Distribution{Univariate, Beta}(Dict(:a=>a, :b=>b))

dims(dist::Distribution{Univariate, Beta}) = ()

vague(::Type{Beta}) = Distribution(Univariate, Beta, a=1.0, b=1.0)

isProper(dist::Distribution{Univariate, Beta}) = (dist.params[:a] > 0.0) && (dist.params[:b] > 0.0)

unsafeMean(dist::Distribution{Univariate, Beta}) = dist.params[:a]/(dist.params[:a] + dist.params[:b])

unsafeLogMean(dist::Distribution{Univariate, Beta}) = digamma(dist.params[:a]) - digamma(dist.params[:a] + dist.params[:b]) # E[log(X)]

unsafeMirroredLogMean(dist::Distribution{Univariate, Beta}) = digamma(dist.params[:b]) - digamma(dist.params[:a] + dist.params[:b]) # E[log(1 - X)]

unsafeVar(dist::Distribution{Univariate, Beta}) = dist.params[:a]*dist.params[:b]/((dist.params[:a] + dist.params[:b])^2*(dist.params[:a] + dist.params[:b] + 1.0))

logPdf(dist::Distribution{Univariate, Beta}, x) = (dist.params[:a]-1)*log(x) + (dist.params[:b]-1)*log(1.0-x) - loggamma(dist.params[:a]) - loggamma(dist.params[:b]) + loggamma(dist.params[:a]+dist.params[:b])

function prod!( x::Distribution{Univariate, Beta},
                y::Distribution{Univariate, Beta},
                z::Distribution{Univariate, Beta}=Distribution(Univariate, Beta, a=1.0, b=1.0))

    z.params[:a] = x.params[:a] + y.params[:a] - 1.0
    z.params[:b] = x.params[:b] + y.params[:b] - 1.0

    return z
end

@symmetrical function prod!(x::Distribution{Univariate, Beta},
                            y::Distribution{Univariate, PointMass},
                            z::Distribution{Univariate, PointMass}=Distribution(Univariate, PointMass, m=0.0))

    (0.0 <= y.params[:m] <= 1.0) || error("PointMass location $(y.params[:m]) should be between 0 and 1")
    z.params[:m] = y.params[:m]

    return z
end

sample(dist::Distribution{Univariate, Beta}) = betainvcdf(dist.params[:a], dist.params[:b], rand())

naturalParams(dist::Distribution{Univariate, Beta}) = [dist.params[:a]-1.0, dist.params[:b]-1.0] # Variant 2 of https://en.wikipedia.org/wiki/Exponential_family

standardDistribution(V::Type{Univariate}, F::Type{Beta}; η::Vector) = Distribution(V, F, a=η[1]+1.0, b=η[2]+1.0)

logNormalizer(::Type{Univariate}, ::Type{Beta}; η::Vector) = loggamma(η[1]+1) + loggamma(η[2]+1) - loggamma(η[1]+η[2]+2)

logPdf(V::Type{Univariate}, F::Type{Beta}, x::Number; η::Vector) = [log(x), log(1-x)]'*η - logNormalizer(V, F; η=η)

# Entropy functional
function differentialEntropy(dist::Distribution{Univariate, Beta})
    labsbeta(dist.params[:a], dist.params[:b]) -
    (dist.params[:a] - 1.0)*digamma(dist.params[:a]) -
    (dist.params[:b] - 1.0)*digamma(dist.params[:b]) +
    (dist.params[:a] + dist.params[:b] - 2.0)*digamma(dist.params[:a] + dist.params[:b])
end

# Average energy functional
function averageEnergy(::Type{Beta}, marg_out::Distribution{Univariate}, marg_a::Distribution{Univariate, PointMass}, marg_b::Distribution{Univariate, PointMass})
    labsbeta(marg_a.params[:m], marg_b.params[:m]) -
    (marg_a.params[:m] - 1.0)*unsafeLogMean(marg_out) -
    (marg_b.params[:m] - 1.0)*unsafeMirroredLogMean(marg_out)
end

# By Stirling's approximation and Monte Carlo summation
function averageEnergy(::Type{Beta}, marg_out::Distribution{Univariate}, marg_a::Distribution{Univariate}, marg_b::Distribution{Univariate})
    unsafeMeanLogMean(marg_a) - unsafeMean(marg_a) + 0.5*(log(2*pi)-unsafeLogMean(marg_a)) +
    unsafeMeanLogMean(marg_b) - unsafeMean(marg_b) + 0.5*(log(2*pi)-unsafeLogMean(marg_b)) -
    (unsafeMean(marg_a)-1)*unsafeLogMean(marg_out) -
    (unsafeMean(marg_b)-1)*unsafeMirroredLogMean(marg_out) -
    sum(loggamma.(sample(marg_a,default_n_samples).+sample(marg_b,default_n_samples)))/default_n_samples
end
