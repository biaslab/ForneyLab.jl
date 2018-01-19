export GaussianMeanPrecision

"""
Description:

    A Gaussian with mean-precision parameterization:

    f(out,m,w) = 𝒩(out|m,w) = (2π)^{-D/2} |w|^{1/2} exp(-1/2 (out - m)' w (out - m))

Interfaces:

    1. out
    2. m (mean)
    3. w (precision)

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
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:m] = self.interfaces[2] = associate!(Interface(self), m)
        self.i[:w] = self.interfaces[3] = associate!(Interface(self), w)

        return self
    end
end

slug(::Type{GaussianMeanPrecision}) = "𝒩"

ProbabilityDistribution(::Type{Univariate}, ::Type{GaussianMeanPrecision}; m=0.0, w=1.0) = ProbabilityDistribution(Univariate, Gaussian, m=m, w=w)
ProbabilityDistribution(::Type{GaussianMeanPrecision}; m::Number=0.0, w::Number=1.0) = ProbabilityDistribution(Univariate, Gaussian, m=m, w=w)
ProbabilityDistribution(::Type{Multivariate}, ::Type{GaussianMeanPrecision}; m=[0.0], w=[1.0].') = ProbabilityDistribution(Multivariate, Gaussian, m=m, w=w)

# Average energy functional
function averageEnergy(::Type{GaussianMeanPrecision}, marg_out::ProbabilityDistribution{Univariate}, marg_mean::ProbabilityDistribution{Univariate}, marg_prec::ProbabilityDistribution{Univariate})
    0.5*log(2*pi) -
    0.5*unsafeLogMean(marg_prec) +
    0.5*unsafeMean(marg_prec)*( unsafeCov(marg_out) + unsafeCov(marg_mean) + (unsafeMean(marg_out) - unsafeMean(marg_mean))^2 )
end

function averageEnergy(::Type{GaussianMeanPrecision}, marg_out::ProbabilityDistribution{Multivariate}, marg_mean::ProbabilityDistribution{Multivariate}, marg_prec::ProbabilityDistribution{MatrixVariate})
    0.5*dims(marg_out)*log(2*pi) -
    0.5*unsafeDetLogMean(marg_prec) +
    0.5*trace( unsafeMean(marg_prec)*(unsafeCov(marg_out) + unsafeCov(marg_mean) + (unsafeMean(marg_out) - unsafeMean(marg_mean))*(unsafeMean(marg_out) - unsafeMean(marg_mean))' ))
end

function averageEnergy(::Type{GaussianMeanPrecision}, marg_out_mean::ProbabilityDistribution{Multivariate, Gaussian}, marg_prec::ProbabilityDistribution{Univariate})
    ensureParameters!(marg_out_mean, (:m, :v))

    V = marg_out_mean.params[:v]
    m = marg_out_mean.params[:m]

    0.5*log(2*pi) -
    0.5*unsafeLogMean(marg_prec) +
    0.5*unsafeMean(marg_prec)*( V[1,1] - 2*V[1,2] + V[2,2] + (m[1] - m[2])^2 )
end

function averageEnergy(::Type{GaussianMeanPrecision}, marg_out_mean::ProbabilityDistribution{Multivariate, Gaussian}, marg_prec::ProbabilityDistribution{MatrixVariate})
    ensureParameters!(marg_out_mean, (:m, :v))

    V = marg_out_mean.params[:v]
    m = marg_out_mean.params[:m]
    d = Int64(dims(marg_out_mean)/2)

    0.5*d*log(2*pi) -
    0.5*unsafeDetLogMean(marg_prec) +
    0.5*trace( unsafeMean(marg_prec)*( V[1:d,1:d] - 2*V[1:d,d+1:end] + V[d+1:end,d+1:end] + (m[1:d] - m[d+1:end])*(m[1:d] - m[d+1:end])' ) )
end