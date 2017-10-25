export GaussianMixture

"""
Description:

    A Gaussian with mean-precision parameterization. For two components:

    f(x, m1, w1, m2, w2, z) = 𝒩(x|m1,w1)^z * 𝒩(x|m1,w1)^(1-z)

Interfaces:

    1. m_1 (mean)
    2. w_1 (precision)
    3. m_2 (mean)
    4. w_2 (precision)
    5. z (switch)
    6. out

Construction:

    GaussianMixture(out, ,m1, m2], [w1, w2], z, id=:some_id)
"""
type GaussianMixture <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function GaussianMixture(out::Variable, m_1::Variable, w_1::Variable, m_2::Variable, w_2::Variable, z::Variable; id=generateId(GaussianMixture))
        self = new(id, Array(Interface, 6), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:m_1] = self.interfaces[1] = associate!(Interface(self), m_1)
        self.i[:w_1] = self.interfaces[2] = associate!(Interface(self), w_1)
        self.i[:m_2] = self.interfaces[3] = associate!(Interface(self), m_2)
        self.i[:w_2] = self.interfaces[4] = associate!(Interface(self), w_2)
        self.i[:z] = self.interfaces[5] = associate!(Interface(self), z)
        self.i[:out] = self.interfaces[6] = associate!(Interface(self), out)

        return self
    end
end

slug(::Type{GaussianMixture}) = "GM"

# Average energy functional
function ForneyLab.averageEnergy(   ::Type{GaussianMixture},
                                    marg_mean_1::Univariate,
                                    marg_prec_1::Univariate, 
                                    marg_mean_2::Univariate,
                                    marg_prec_2::Univariate, 
                                    marg_switch::Univariate{Bernoulli},
                                    marg_out::Univariate)

    unsafeMean(marg_switch)*averageEnergy(GaussianMeanPrecision, marg_mean_1, marg_prec_1, marg_out) +
    (1.0 - unsafeMean(marg_switch))*averageEnergy(GaussianMeanPrecision, marg_mean_2, marg_prec_2, marg_out)
end