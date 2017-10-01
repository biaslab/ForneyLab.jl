export GaussianMeanVariance

"""
Description:

    A Gaussian with mean-variance parameterization:

    f(x,m,v) = 𝒩(x|m,v)

Interfaces:

    1. mean
    2. variance
    3. out

Construction:

    GaussianMeanVariance(out, mean, variance, id=:some_id)
"""
type GaussianMeanVariance <: Gaussian
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function GaussianMeanVariance(out::Variable, mean::Variable, variance::Variable; id=generateId(Gaussian))
        self = new(id, Array(Interface, 3), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:mean] = self.interfaces[1] = associate!(Interface(self), mean)
        self.i[:variance] = self.interfaces[2] = associate!(Interface(self), variance)
        self.i[:out] = self.interfaces[3] = associate!(Interface(self), out)

        return self
    end
end

slug(::Type{GaussianMeanVariance}) = "𝒩"

# Average energy functional
function averageEnergy(::Type{GaussianMeanVariance}, marg_mean::ProbabilityDistribution, marg_var::ProbabilityDistribution, marg_out::ProbabilityDistribution)
    0.5*log(2*pi) +
    0.5*unsafeLogMean(marg_var) +
    0.5*unsafeInverseMean(marg_var)*( unsafeCov(marg_out) + unsafeCov(marg_mean) + (unsafeMean(marg_out) - unsafeMean(marg_mean))^2 )
end
