export Gamma

"""
Description:

    A gamma node with shape-rate parameterization:

    f(x,a,b) = Gam(x|a,b)

Interfaces:

    1. a (shape)
    2. b (rate)
    3. out

Construction:

    Gamma(out, a, b, id=:some_id)
"""
type Gamma <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function Gamma(out::Variable, a::Variable, b::Variable; id=generateId(Gamma))
        self = new(id, Array(Interface, 3), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:a] = self.interfaces[1] = associate!(Interface(self), a)
        self.i[:b] = self.interfaces[2] = associate!(Interface(self), b)
        self.i[:out] = self.interfaces[3] = associate!(Interface(self), out)

        return self
    end
end

slug(::Type{Gamma}) = "Gam"

Univariate(::Type{Gamma}; a=1.0, b=1.0) = Univariate{Gamma}(Dict(a=>1.0, b=>1.0))

vague(::Type{Univariate{Gamma}}) = Univariate(Gamma, a=1.0, b=tiny) # Flat prior leads to more stable behaviour than Jeffrey's prior

unsafeMean(dist::Univariate{Gamma}) = dist.params[:a]/dist.params[:b] # unsafe mean

unsafeLogMean(dist::Univariate{Gamma}) = digamma(dist.params[:a]) - log(dist.params[:b])

unsafeVar(dist::Univariate{Gamma}) = dist.params[:a]/dist.params[:b]^2 # unsafe variance

isProper(dist::Univariate{Gamma}) = (dist.params[:a] >= tiny) && (dist.params[:b] >= tiny)

function prod!( x::Univariate{Gamma},
                y::Univariate{Gamma},
                z::Univariate{Gamma}=Univariate(Gamma, a=0.0, b=0.0))

    z.params[:a] = x.params[:a] + y.params[:a] - 1.0
    z.params[:b] = x.params[:b] + y.params[:b]

    return z
end

# Entropy functional
function differentialEntropy(dist::Univariate{Gamma})
    lgamma(dist.params[:a]) -
    (dist.params[:a] - 1.0)*digamma(dist.params[:a]) -
    log(dist.params[:b]) +
    dist.params[:a]
end

# Average energy functional
function averageEnergy(::Type{Gamma}, marg_a::Univariate{PointMass}, marg_b::Univariate, marg_out::Univariate)
    lgamma(marg_a.params[:m]) -
    marg_a.params[:m]*unsafeLogMean(marg_b) -
    (marg_a.params[:m] - 1.0)*unsafeLogMean(marg_out) +
    unsafeMean(marg_b)*unsafeMean(marg_out)
end