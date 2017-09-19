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

    function GaussianMeanVariance(out::Variable, mean::Variable, variance::Variable; id=generateId(GaussianMeanVariance))
        self = new(id, Array(Interface, 3), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:mean] = self.interfaces[1] = associate!(Interface(self), mean)
        self.i[:variance] = self.interfaces[2] = associate!(Interface(self), variance)
        self.i[:out] = self.interfaces[3] = associate!(Interface(self), out)

        return self
    end
end

slug(::Type{GaussianMeanVariance}) = "𝒩"

ProbabilityDistribution(::Type{GaussianMeanVariance}) = ProbabilityDistribution(GaussianMeanVariance, m=0.0, v=1.0)

# TODO: make more efficient by introducing alternative Gaussian parameterizations
function prod!( x::ProbabilityDistribution{GaussianMeanVariance},
                y::ProbabilityDistribution{GaussianMeanVariance},
                z::ProbabilityDistribution{GaussianMeanVariance}=ProbabilityDistribution(GaussianMeanVariance))

    # Multiplication of 2 Gaussian PDFs: p(z) = p(x) * p(y)
    w_x = 1/x.params[:v]
    xi_x = w_x*x.params[:m]
    w_y = 1/y.params[:v]
    xi_y = w_y*y.params[:m]
    w_z = w_x + w_y
    xi_z = xi_x + xi_y

    z.params[:m] = xi_z/w_z
    z.params[:v] = 1/w_z

    return z
end