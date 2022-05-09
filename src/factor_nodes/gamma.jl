export Gamma

"""
Description:

    A gamma node with shape-rate parameterization:

    f(out,a,b) = Gam(out|a,b) = 1/Γ(a) b^a out^{a - 1} exp(-b out)

Interfaces:

    1. out
    2. a (shape)
    3. b (rate)

Construction:

    Gamma(out, a, b, id=:some_id)
"""
mutable struct Gamma <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function Gamma(out, a, b; id=generateId(Gamma))
        @ensureVariables(out, a, b)
        self = new(id, Array{Interface}(undef, 3), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:a] = self.interfaces[2] = associate!(Interface(self), a)
        self.i[:b] = self.interfaces[3] = associate!(Interface(self), b)

        return self
    end
end

slug(::Type{Gamma}) = "Gam"

format(dist::Distribution{Univariate, Gamma}) = "$(slug(Gamma))(a=$(format(dist.params[:a])), b=$(format(dist.params[:b])))"

Distribution(::Type{Univariate}, ::Type{Gamma}; a=1.0, b=1.0) = Distribution{Univariate, Gamma}(Dict(:a=>a, :b=>b))
Distribution(::Type{Gamma}; a=1.0, b=1.0) = Distribution{Univariate, Gamma}(Dict(:a=>a, :b=>b))

dims(dist::Distribution{Univariate, Gamma}) = ()

vague(::Type{Gamma}) = Distribution(Univariate, Gamma, a=1.0, b=tiny) # Flat prior leads to more stable behaviour than Jeffrey's prior

unsafeMean(dist::Distribution{Univariate, Gamma}) = dist.params[:a]/dist.params[:b] # unsafe mean

unsafeLogMean(dist::Distribution{Univariate, Gamma}) = digamma(dist.params[:a]) - log(dist.params[:b])

# https://stats.stackexchange.com/questions/457357/what-is-the-expected-value-of-x-logx-of-the-gamma-distribution
unsafeMeanLogMean(dist::Distribution{Univariate, Gamma}) =  (gamma(dist.params[:a]+1)/(gamma(dist.params[:a])*dist.params[:b])) * (digamma(dist.params[:a]+1) - log(dist.params[:b]))

unsafeVar(dist::Distribution{Univariate, Gamma}) = dist.params[:a]/dist.params[:b]^2 # unsafe variance

logPdf(dist::Distribution{Univariate, Gamma}, x) = dist.params[:a]*log(dist.params[:b]) - labsgamma(dist.params[:a]) + (dist.params[:a]-1)*log(x) - dist.params[:b]*x

isProper(dist::Distribution{Univariate, Gamma}) = (dist.params[:a] >= tiny) && (dist.params[:b] >= tiny)

function prod!( x::Distribution{Univariate, Gamma},
                y::Distribution{Univariate, Gamma},
                z::Distribution{Univariate, Gamma}=Distribution(Univariate, Gamma, a=0.0, b=0.0))

    z.params[:a] = x.params[:a] + y.params[:a] - 1.0
    z.params[:b] = x.params[:b] + y.params[:b]

    return z
end

@symmetrical function prod!(x::Distribution{Univariate, Gamma},
                            y::Distribution{Univariate, PointMass},
                            z::Distribution{Univariate, PointMass}=Distribution(Univariate, PointMass, m=0.0))

    (y.params[:m] > 0.0) || error("PointMass location $(y.params[:m]) should be positive")
    z.params[:m] = y.params[:m]

    return z
end

sample(dist::Distribution{Univariate, Gamma}) = gammainvcdf(dist.params[:a], 1/dist.params[:b], rand())

naturalParams(dist::Distribution{Univariate, Gamma}) = [dist.params[:a]-1.0, -dist.params[:b]]

standardDistribution(V::Type{Univariate}, F::Type{Gamma}; η::Vector) = Distribution(V, F, a=η[1]+1.0, b=-η[2])

logNormalizer(::Type{Univariate}, ::Type{Gamma}; η::Vector) = loggamma(η[1]+1.0) - (η[1]+1.0)*log(-η[2])

logPdf(V::Type{Univariate}, F::Type{Gamma}, x::Number; η::Vector) =[log(x), x]'*η - logNormalizer(V, F, η=η)

# Entropy functional
function differentialEntropy(dist::Distribution{Univariate, Gamma})
    labsgamma(dist.params[:a]) -
    (dist.params[:a] - 1.0)*digamma(dist.params[:a]) -
    log(dist.params[:b]) +
    dist.params[:a]
end

# Average energy functional
function averageEnergy(::Type{Gamma}, marg_out::Distribution{Univariate}, marg_a::Distribution{Univariate, PointMass}, marg_b::Distribution{Univariate})
    labsgamma(marg_a.params[:m]) -
    marg_a.params[:m]*unsafeLogMean(marg_b) -
    (marg_a.params[:m] - 1.0)*unsafeLogMean(marg_out) +
    unsafeMean(marg_b)*unsafeMean(marg_out)
end

# By Stirling's approximation
function averageEnergy(::Type{Gamma}, marg_out::Distribution{Univariate}, marg_a::Distribution{Univariate}, marg_b::Distribution{Univariate})
    unsafeMeanLogMean(marg_a) - unsafeMean(marg_a) + 0.5*(log(2*pi)-unsafeLogMean(marg_a)) -
    unsafeMean(marg_a)*unsafeLogMean(marg_b) -
    (unsafeMean(marg_a) - 1.0)*unsafeLogMean(marg_out) +
    unsafeMean(marg_b)*unsafeMean(marg_out)
end
