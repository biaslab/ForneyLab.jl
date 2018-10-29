export GaussianMixture

"""
Description:

    A Gaussian mixture with mean-precision parameterization:

    f(out, z, m1, w1, m2, w2, ...) = 𝒩(out|m1, w1)^z_1 * 𝒩(out|m2, w2)^z_2 * ...

Interfaces:

    1. out
    2. z (switch)
    3. m1 (mean)
    4. w1 (precision)
    5. m2 (mean)
    6. w2 (precision)
    ...

Construction:

    GaussianMixture(out, z, m1, w1, m2, w2, ..., id=:some_id)
"""
mutable struct GaussianMixture <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function GaussianMixture(out, z, args::Vararg; id=generateId(GaussianMixture))
        n_args = length(args)
        @ensureVariables(out, z)
        for i=1:n_args
            @ensureVariables(args[i])
        end
        iseven(n_args) || error("Number of mixture arguments should be even")
        self = new(id, Array{Interface}(undef, length(args) + 2), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:z] = self.interfaces[2] = associate!(Interface(self), z)
        n_factors = Int64(n_args/2)
        for k = 1:n_factors
            self.i[:m*k] = self.interfaces[2*k + 1] = associate!(Interface(self), args[2*k - 1])
            self.i[:w*k] = self.interfaces[2*k + 2] = associate!(Interface(self), args[2*k])
        end

        return self
    end
end

slug(::Type{GaussianMixture}) = "GM"

# Average energy functional
function ForneyLab.averageEnergy(   ::Type{GaussianMixture},
                                    dist_out::ProbabilityDistribution,
                                    dist_switch::ProbabilityDistribution,
                                    dist_factors::Vararg{ProbabilityDistribution})

    dist_means = collect(dist_factors[1:2:end])
    dist_precs = collect(dist_factors[2:2:end])
    n_factors = length(dist_means)
    z_bar = unsafeMeanVector(dist_switch)

    U = 0.0
    for k = 1:n_factors
        U += z_bar[k]*averageEnergy(GaussianMeanPrecision, dist_out, dist_means[k], dist_precs[k])
    end

    return U
end