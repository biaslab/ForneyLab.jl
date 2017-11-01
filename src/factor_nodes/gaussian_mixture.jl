export GaussianMixture

"""
Description:

    A Gaussian with mean-precision parameterization. For two components:

    f(out, m1, w1, m2, w2, z) = 𝒩(out|m1,w1)^z * 𝒩(out|m1,w1)^(1-z)

Interfaces:

    1. out
    2. z (switch)
    3. m1 (mean)
    4. w1 (precision)
    5. m2 (mean)
    6. w2 (precision)

Construction:

    GaussianMixture(out, z, m1, m2, w1, w2, id=:some_id)
"""
type GaussianMixture <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function GaussianMixture(out::Variable, z::Variable, m1::Variable, w1::Variable, m2::Variable, w2::Variable; id=generateId(GaussianMixture))
        self = new(id, Array(Interface, 6), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:z] = self.interfaces[2] = associate!(Interface(self), z)
        self.i[:m1] = self.interfaces[3] = associate!(Interface(self), m1)
        self.i[:w1] = self.interfaces[4] = associate!(Interface(self), w1)
        self.i[:m2] = self.interfaces[5] = associate!(Interface(self), m2)
        self.i[:w2] = self.interfaces[6] = associate!(Interface(self), w2)

        return self
    end
end

slug(::Type{GaussianMixture}) = "GM"

# Average energy functional
function ForneyLab.averageEnergy(   ::Type{GaussianMixture},
                                    marg_out::ProbabilityDistribution{Univariate},
                                    marg_switch::ProbabilityDistribution{Univariate, Bernoulli},
                                    marg_mean_1::ProbabilityDistribution{Univariate},
                                    marg_prec_1::ProbabilityDistribution{Univariate}, 
                                    marg_mean_2::ProbabilityDistribution{Univariate},
                                    marg_prec_2::ProbabilityDistribution{Univariate})

    unsafeMean(marg_switch)*averageEnergy(GaussianMeanPrecision, marg_out, marg_mean_1, marg_prec_1) +
    (1.0 - unsafeMean(marg_switch))*averageEnergy(GaussianMeanPrecision, marg_out, marg_mean_2, marg_prec_2)
end