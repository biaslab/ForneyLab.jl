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
    ...

Construction:

    GaussianMixture(out, z, m1, m2, w1, w2, ..., id=:some_id)
"""
type GaussianMixture <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function GaussianMixture(out::Variable, z::Variable, args::Vararg{Variable}; id=generateId(GaussianMixture))
        n = length(args)
        iseven(n) || error("Number of mixture variables should be even")
        self = new(id, Array(Interface, length(args) + 2), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:z] = self.interfaces[2] = associate!(Interface(self), z)
        for c = 1:Int64(n/2)
            self.i[:m*c] = self.interfaces[2*c + 1] = associate!(Interface(self), args[2*c - 1])
            self.i[:w*c] = self.interfaces[2*c + 2] = associate!(Interface(self), args[2*c])
        end

        return self
    end
end

slug(::Type{GaussianMixture}) = "GM"

# Average energy functional
function ForneyLab.averageEnergy(   ::Type{GaussianMixture},
                                    marg_out::ProbabilityDistribution,
                                    marg_switch::ProbabilityDistribution{Univariate, Bernoulli},
                                    marg_mean_1::ProbabilityDistribution,
                                    marg_prec_1::ProbabilityDistribution, 
                                    marg_mean_2::ProbabilityDistribution,
                                    marg_prec_2::ProbabilityDistribution)

    unsafeMean(marg_switch)*averageEnergy(GaussianMeanPrecision, marg_out, marg_mean_1, marg_prec_1) +
    (1.0 - unsafeMean(marg_switch))*averageEnergy(GaussianMeanPrecision, marg_out, marg_mean_2, marg_prec_2)
end