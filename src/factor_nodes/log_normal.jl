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

format(dist::Distribution{Univariate, LogNormal}) = "$(slug(LogNormal))(m=$(format(dist.params[:m])), s=$(format(dist.params[:s])))"

Distribution(::Type{Univariate}, ::Type{LogNormal}; m::Float64=1.0, s::Float64=1.0) = Distribution{Univariate, LogNormal}(Dict(:m=>m, :s=>s))
Distribution(::Type{LogNormal}; m::Float64=1.0, s::Float64=1.0) = Distribution{Univariate, LogNormal}(Dict(:m=>m, :s=>s))

dims(dist::Distribution{Univariate, LogNormal}) = ()

vague(::Type{LogNormal}) = Distribution(Univariate, LogNormal, m=1.0, s=huge)

unsafeMean(dist::Distribution{Univariate, LogNormal}) = exp(dist.params[:m] + 0.5*dist.params[:s])
unsafeLogMean(dist::Distribution{Univariate, LogNormal}) = dist.params[:m]

unsafeVar(dist::Distribution{Univariate, LogNormal}) = (exp(dist.params[:s]) - 1.0)*exp(2.0*dist.params[:m] + dist.params[:s])
unsafeLogVar(dist::Distribution{Univariate, LogNormal}) = dist.params[:s]

unsafeCov(dist::Distribution{Univariate, LogNormal}) = unsafeVar(dist)
unsafeLogCov(dist::Distribution{Univariate, LogNormal}) = dist.params[:s]

logPdf(dist::Distribution{Univariate, LogNormal},x) = -0.5*(log(2pi)+log(dist.params[:s])) -log(x) -0.5*(log(x)-dist.params[:m])^2/dist.params[:s]
isProper(dist::Distribution{Univariate, LogNormal}) = (dist.params[:s] > 0.0)

sample(dist::Distribution{Univariate, LogNormal}) = exp(dist.params[:m]+sqrt(dist.params[:s])*randn())

naturalParams(dist::Distribution{Univariate, LogNormal}) = [dist.params[:m]/dist.params[:s], -0.5/dist.params[:s]]

standardDistribution(V::Type{Univariate}, F::Type{LogNormal}; η::Vector) = Distribution(V, F, m=-(0.5/η[2])*η[1], s=-0.5/η[2])

logNormalizer(::Type{Univariate}, ::Type{LogNormal}; η::Vector) = -η[1]^2/(4*η[2]) - 0.5*log(-2*η[2])

logPdf(V::Type{Univariate}, F::Type{LogNormal}, x::Number; η::Vector) = -0.5*log(2pi) - log(x) + [log(x), log(x)^2]'*η - logNormalizer(V, F, η=η)

"""
Gamma approximation to the log-normal distribution using Laplace's method
"""
laplace(::Type{Gamma}, dist::Distribution{Univariate, LogNormal}) = Distribution(Univariate, Gamma, a=1/dist.params[:s], b=1/dist.params[:s]*exp(-dist.params[:m]))

@symmetrical function prod!(x::Distribution{Univariate, LogNormal},
                            y::Distribution{Univariate, Gamma},
                            z::Distribution{Univariate, Gamma}=Distribution(Univariate, Gamma, a=1.0, b=1.0))

    x_approx = laplace(Gamma, x)
    z.params[:a] = x_approx.params[:a] + y.params[:a] - 1.0
    z.params[:b] = x_approx.params[:b] + y.params[:b]

    return z
end

@symmetrical function prod!(x::Distribution{Univariate, LogNormal},
                            y::Distribution{Univariate, PointMass},
                            z::Distribution{Univariate, PointMass}=Distribution(Univariate, PointMass, m=0.0))

    (y.params[:m] > 0.0) || error("PointMass location $(y.params[:m]) should be positive")
    z.params[:m] = y.params[:m]

    return z
end

# Entropy functional
function differentialEntropy(dist::Distribution{Univariate, LogNormal})
    0.5*log(dist.params[:s]) +
    dist.params[:m] + 0.5 +
    0.5*log(2*pi)
end

# Average energy functional
function averageEnergy(::Type{LogNormal}, marg_out::Distribution{Univariate}, marg_m::Distribution{Univariate}, marg_s::Distribution{Univariate})
    unsafeLogMean(marg_out) +
    0.5*log(2*pi) +
    0.5*unsafeLogMean(marg_s) +
    0.5*unsafeInverseMean(marg_s)*( unsafeCov(marg_m) + unsafeLogCov(marg_out) + (unsafeLogMean(marg_out) - unsafeMean(marg_m))^2 )
end
