export TransitionMixture

"""
Description:

    A mixture of transition matrices:

    f(out, z, A1, A2, ...) = Dir(out|A1)^z_1 * Dir(out|A2)^z_2 * ...

Interfaces:

    1. out
    2. z (switch)
    3. A1 (transition matrix)
    4. A2 (transition matrix)
    ...

Construction:

    TransitionMixture(out, z, A1, A2, ..., id=:some_id)
"""
mutable struct TransitionMixture <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function TransitionMixture(out, z, args::Vararg; id=generateId(TransitionMixture))
        n_factors = length(args)
        @ensureVariables(out, z)
        for i=1:n_factors
            @ensureVariables(args[i])
        end
        self = new(id, Array{Interface}(undef, n_factors + 2), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:z] = self.interfaces[2] = associate!(Interface(self), z)
        for k = 1:n_factors
            self.i[:a*k] = self.interfaces[k + 2] = associate!(Interface(self), args[k])
        end

        return self
    end
end

slug(::Type{TransitionMixture}) = "TM"

# Average energy functional
function ForneyLab.averageEnergy(::Type{TransitionMixture},
                                 dist_out::ProbabilityDistribution,
                                 dist_switch::ProbabilityDistribution,
                                 dist_factors::Vararg{ProbabilityDistribution{MatrixVariate, PointMass}})

    n_factors = length(dist_factors)
    z_bar = unsafeMeanVector(dist_switch)
    U = 0.0
    for k = 1:n_factors
        A_k = clamp.(dist_factors[k].params[:m], tiny, 1-tiny) # Soften given transition functions
        U += z_bar[k]*averageEnergy(Dirichlet, dist_out, ProbabilityDistribution(MatrixVariate, PointMass, m=A_k))
    end

    return U
end