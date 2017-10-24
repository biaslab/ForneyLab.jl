export GaussianMeanPrecision

"""
Description:

    A Gaussian with mean-precision parameterization:

    f(x,m,w) = 𝒩(x|m,w)

Interfaces:

    1. m (mean)
    2. w (precision)
    3. out

Construction:

    GaussianMeanPrecision(out, m, w, id=:some_id)
"""
type GaussianMeanPrecision <: Gaussian
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function GaussianMeanPrecision(out::Variable, m::Variable, w::Variable; id=generateId(Gaussian))
        self = new(id, Array(Interface, 3), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:m] = self.interfaces[1] = associate!(Interface(self), m)
        self.i[:w] = self.interfaces[2] = associate!(Interface(self), w)
        self.i[:out] = self.interfaces[3] = associate!(Interface(self), out)

        return self
    end
end

slug(::Type{GaussianMeanPrecision}) = "𝒩"

# Average energy functional
function averageEnergy(::Type{GaussianMeanPrecision}, marg_mean::Univariate, marg_prec::Univariate, marg_out::Univariate)
    0.5*log(2*pi) -
    0.5*unsafeLogMean(marg_prec) +
    0.5*unsafeMean(marg_prec)*( unsafeCov(marg_out) + unsafeCov(marg_mean) + (unsafeMean(marg_out) - unsafeMean(marg_mean))^2 )
end

function averageEnergy{dims}(::Type{GaussianMeanPrecision}, marg_mean::Multivariate{dims}, marg_prec::MatrixVariate{dims, dims}, marg_out::Multivariate{dims})
    0.5*dims*log(2*pi) -
    0.5*unsafeDetLogMean(marg_prec) +
    0.5*trace( unsafeMean(marg_prec)*(unsafeCov(marg_out) + unsafeCov(marg_mean) + (unsafeMean(marg_out) - unsafeMean(marg_mean))*(unsafeMean(marg_out) - unsafeMean(marg_mean))' ))
end