export ChanceConstraint

"""
Description:

    Chance constraint on the marginal q of the connected variable x, as
    ϵ ⩽ 1 - ∫_G q(x) dx, where G indicates a region. In other words,
    the probability mass of q is not allowed to overflow the region G
    by more than ϵ. Implementation according to (van de Laar et al.
    "Chance-Constrained Active Inference", MIT Neural Computation, 2021).
    
Interfaces:

    1. out

Construction:

    ChanceConstraint(out; G=(min,max), epsilon=epsilon, id=:my_node)
"""
mutable struct ChanceConstraint <: SoftFactor
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol, Interface}

    G::Tuple # Range for G (1-D only)
    epsilon::Float64 # Threshold
    atol::Union{Float64, Nothing} # Tolerance for corrected belief

    function ChanceConstraint(out; G::Tuple, epsilon::Float64, atol=nothing, id=ForneyLab.generateId(ChanceConstraint))
        @ensureVariables(out)
        self = new(id, Vector{Interface}(undef, 1), Dict{Symbol,Interface}(), G, epsilon, atol)
        ForneyLab.addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)

        return self
    end
end    

slug(::Type{ChanceConstraint}) = "P"

# A breaker message is required if interface is partnered with an expectation constraint
requiresBreaker(interface::Interface, partner_interface::Interface, partner_node::ChanceConstraint) = true

breakerParameters(interface::Interface, partner_interface::Interface, partner_node::ChanceConstraint) = (Message{GaussianMeanVariance, Univariate}, ()) # Univariate only

# Constraints do not contribute to average energy
averageEnergy(::Type{ChanceConstraint}, marg_out::ProbabilityDistribution) = 0.0
