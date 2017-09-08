export GaussianMeanPrecision

"""
Description:

    A Gaussian with mean-precision parameterization:

    f(x,m,w) = 𝒩(x|m,w)

Interfaces:

    1. mean
    2. precision
    3. out

Construction:

    GaussianMeanPrecision(out, mean, precision, id=:some_id)
"""
type GaussianMeanPrecision <: Gaussian
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function GaussianMeanPrecision(out::Variable, mean::Variable, precision::Variable; id=generateId(Gaussian))
        self = new(id, Array(Interface, 3), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:mean] = self.interfaces[1] = associate!(Interface(self), mean)
        self.i[:precision] = self.interfaces[2] = associate!(Interface(self), precision)
        self.i[:out] = self.interfaces[3] = associate!(Interface(self), out)

        return self
    end
end

slug(::Type{GaussianMeanPrecision}) = "𝒩"