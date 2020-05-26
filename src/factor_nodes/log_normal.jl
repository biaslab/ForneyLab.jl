export LogNormal, momentMatching, logMomentMatching

"""
Description:

    A log-normal node with location-scale parameterization:

    f(out,m,s) = logN(out|m, s) = 1/out (2π s)^{-1/2} exp(-1/(2s) (log(out) - m)^2))

Interfaces:

    1. out
    2. m (location)
    3. s (squared scale)

Construction:

    LogNormal(out, m, s, id=:some_id)
"""
mutable struct LogNormal <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function LogNormal(out, m, s; id=generateId(LogNormal))
        @ensureVariables(out, m, s)
        self = new(id, Array{Interface}(undef, 3), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:m] = self.interfaces[2] = associate!(Interface(self), m)
        self.i[:s] = self.interfaces[3] = associate!(Interface(self), s)

        return self
    end
end

slug(::Type{LogNormal}) = "log𝒩"

format(dist::ProbabilityDistribution{Univariate, LogNormal}) = "$(slug(LogNormal))(m=$(format(dist.params[:m])), s=$(format(dist.params[:s])))"

ProbabilityDistribution(::Type{Univariate}, ::Type{LogNormal}; m::Float64=1.0, s::Float64=1.0) = ProbabilityDistribution{Univariate, LogNormal}(Dict(:m=>m, :s=>s))
ProbabilityDistribution(::Type{LogNormal}; m::Float64=1.0, s::Float64=1.0) = ProbabilityDistribution{Univariate, LogNormal}(Dict(:m=>m, :s=>s))

dims(dist::ProbabilityDistribution{Univariate, LogNormal}) = 1

vague(::Type{LogNormal}) = ProbabilityDistribution(Univariate, LogNormal, m=1.0, s=huge)

unsafeMean(dist::ProbabilityDistribution{Univariate, LogNormal}) = exp(dist.params[:m] + 0.5*dist.params[:s])
unsafeLogMean(dist::ProbabilityDistribution{Univariate, LogNormal}) = dist.params[:m]

unsafeVar(dist::ProbabilityDistribution{Univariate, LogNormal}) = (exp(dist.params[:s]) - 1.0)*exp(2.0*dist.params[:m] + dist.params[:s])
unsafeLogVar(dist::ProbabilityDistribution{Univariate, LogNormal}) = dist.params[:s]

unsafeCov(dist::ProbabilityDistribution{Univariate, LogNormal}) = unsafeVar(dist)
unsafeLogCov(dist::ProbabilityDistribution{Univariate, LogNormal}) = dist.params[:s]

logPdf(dist::ProbabilityDistribution{Univariate, LogNormal},x) = -0.5*(log(2pi)+log(dist.params[:s])) -log(x) -0.5*(log(x)-dist.params[:m])^2/dist.params[:s]
isProper(dist::ProbabilityDistribution{Univariate, LogNormal}) = (dist.params[:s] > 0.0)

function sample(dist::ProbabilityDistribution{Univariate, LogNormal})
    return exp(dist.params[:m]+sqrt(dist.params[:s])*randn())
end

"""
Gamma approximation to the log-normal distribution using Laplace's method
"""
laplace(::Type{Gamma}, dist::ProbabilityDistribution{Univariate, LogNormal}) = ProbabilityDistribution(Univariate, Gamma, a=1/dist.params[:s], b=1/dist.params[:s]*exp(-dist.params[:m]))

@symmetrical function prod!(x::ProbabilityDistribution{Univariate, LogNormal},
                            y::ProbabilityDistribution{Univariate, Gamma},
                            z::ProbabilityDistribution{Univariate, Gamma}=ProbabilityDistribution(Univariate, Gamma, a=1.0, b=1.0))

    x_approx = laplace(Gamma, x)
    z.params[:a] = x_approx.params[:a] + y.params[:a] - 1.0
    z.params[:b] = x_approx.params[:b] + y.params[:b]

    return z
end

@symmetrical function prod!(x::ProbabilityDistribution{Univariate, LogNormal},
                            y::ProbabilityDistribution{Univariate, PointMass},
                            z::ProbabilityDistribution{Univariate, PointMass}=ProbabilityDistribution(Univariate, PointMass, m=0.0))

    (y.params[:m] > 0.0) || error("PointMass location $(y.params[:m]) should be positive")
    z.params[:m] = y.params[:m]

    return z
end

# Entropy functional
function differentialEntropy(dist::ProbabilityDistribution{Univariate, LogNormal})
    0.5*log(dist.params[:s]) +
    dist.params[:m] + 0.5 +
    0.5*log(2*pi)
end

# Average energy functional
function averageEnergy(::Type{LogNormal}, marg_out::ProbabilityDistribution{Univariate}, marg_m::ProbabilityDistribution{Univariate}, marg_s::ProbabilityDistribution{Univariate})
    unsafeLogMean(marg_out) +
    0.5*log(2*pi) +
    0.5*unsafeLogMean(marg_s) +
    0.5*unsafeInverseMean(marg_s)*( unsafeCov(marg_m) + unsafeLogCov(marg_out) + (unsafeLogMean(marg_out) - unsafeMean(marg_m))^2 )
end
