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
type GaussianMeanVariance <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function GaussianMeanVariance(out::Variable, mean::Variable, variance::Variable; id=generateId(GaussianMeanVariance))
        self = new(id, Array(Interface, 3), Dict{Symbol,Interface}())
        self.i[:mean] = self.interfaces[1] = associate!(Interface(self), mean)
        self.i[:variance] = self.interfaces[2] = associate!(Interface(self), variance)
        self.i[:out] = self.interfaces[3] = associate!(Interface(self), out)

        addNode!(currentGraph(), self)

        return out
    end
end

# TODO: rules id dictionary