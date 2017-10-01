export Gamma

"""
Description:

    A gamma node with shape-rate parameterization:

    f(x,a,b) = Gam(x|a,b)

Interfaces:

    1. shape
    2. rate
    3. out

Construction:

    Gamma(out, shape, rate, id=:some_id)
"""
type Gamma <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function Gamma(out::Variable, shape::Variable, rate::Variable; id=generateId(Gamma))
        self = new(id, Array(Interface, 3), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:shape] = self.interfaces[1] = associate!(Interface(self), shape)
        self.i[:rate] = self.interfaces[2] = associate!(Interface(self), rate)
        self.i[:out] = self.interfaces[3] = associate!(Interface(self), out)

        return self
    end
end

slug(::Type{Gamma}) = "Gam"

ProbabilityDistribution(::Type{Gamma}) = ProbabilityDistribution(Gamma, a=1.0, b=1.0)

vague(::Type{ProbabilityDistribution{Gamma}}) = ProbabilityDistribution(Gamma, a=tiny, b=tiny)

unsafeMean(dist::ProbabilityDistribution{Gamma}) = dist.params[:a]/dist.params[:b] # unsafe mean

unsafeLogMean(dist::ProbabilityDistribution{Gamma}) = digamma(dist.params[:a]) - log(dist.params[:b])

unsafeVar(dist::ProbabilityDistribution{Gamma}) = dist.params[:a]/dist.params[:b]^2 # unsafe variance

isProper(dist::ProbabilityDistribution{Gamma}) = (dist.params[:a] >= tiny) && (dist.params[:b] >= tiny)

function prod!( x::ProbabilityDistribution{Gamma},
                y::ProbabilityDistribution{Gamma},
                z::ProbabilityDistribution{Gamma}=ProbabilityDistribution(Gamma, a=0.0, b=0.0))

    z.params[:a] = x.params[:a] + y.params[:a] - 1.0
    z.params[:b] = x.params[:b] + y.params[:b]

    return z
end

# Entropy functional
function differentialEntropy(dist::ProbabilityDistribution{Gamma})
    lgamma(dist.params[:a]) -
    (dist.params[:a] - 1.0)*digamma(dist.params[:a]) -
    log(dist.params[:b]) +
    dist.params[:a]
end

# Average energy functional
function averageEnergy(::Type{Gamma}, marg_a::ProbabilityDistribution{PointMass}, marg_b::ProbabilityDistribution, marg_out::ProbabilityDistribution)
    lgamma(marg_a.params[:m]) -
    marg_a.params[:m]*unsafeLogMean(marg_b) -
    (marg_a.params[:m] - 1.0)*unsafeLogMean(marg_out) +
    unsafeMean(marg_b)*unsafeMean(marg_out)
end