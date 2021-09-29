export TransitionMixture

"""
Description:

    A mixture of discrete transitions:

    f(out, in1, z, A1, A2, ...) = Cat(out|A1*in1)^z_1 * Cat(out|A2*in1)^z_2 * ...

Interfaces:

    1. out
    2. in1
    3. z (switch)
    4. A1 (transition matrix)
    5. A2 (transition matrix)
    ...

Construction:

    TransitionMixture(out, in1, z, A1, A2, ..., id=:some_id)
"""
mutable struct TransitionMixture <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function TransitionMixture(out, in1, z, args::Vararg; id=generateId(TransitionMixture))
        n_factors = length(args)
        @ensureVariables(out, in1, z)
        for i=1:n_factors
            @ensureVariables(args[i])
        end
        self = new(id, Array{Interface}(undef, n_factors + 3), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:in1] = self.interfaces[2] = associate!(Interface(self), in1)
        self.i[:z] = self.interfaces[3] = associate!(Interface(self), z)
        for k = 1:n_factors
            self.i[:a*k] = self.interfaces[k + 3] = associate!(Interface(self), args[k])
        end

        return self
    end
end

slug(::Type{TransitionMixture}) = "TM"

# Average energy functional
function averageEnergy(::Type{TransitionMixture},
                       dist_out_in1_switch::ProbabilityDistribution{Multivariate, Contingency},
                       dist_factors::Vararg{ProbabilityDistribution})

    n_factors = length(dist_factors)
    U = 0.0
    for k = 1:n_factors
        U += -tr(dist_out_in1_switch.params[:p][k]'*unsafeLogMean(dist_factors[k]))
    end

    return U
end