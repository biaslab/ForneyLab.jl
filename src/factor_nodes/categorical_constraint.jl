export CategoricalConstraint

"""
Description:

    Constraints the marginal of the connected variable to a Categorical belief.
    
Interfaces:

    1. out

Construction:

    CategoricalConstraint(out; p=p, id=:my_node)
"""
mutable struct CategoricalConstraint <: SoftFactor # TODO: how to handle free energy evaluation for a point-mass constrained variable?
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol, Interface}

    p::Vector{Float64}

    function CategoricalConstraint(out; p::Vector{Float64},  id=ForneyLab.generateId(CategoricalConstraint))
        @ensureVariables(out)
        self = new(id, Vector{Interface}(undef, 1), Dict{Symbol,Interface}(), p)
        ForneyLab.addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)

        return self
    end
end    

slug(::Type{CategoricalConstraint}) = "Cat_q"

# A breaker message is required if interface is partnered with a point-mass constraint
requiresBreaker(interface::Interface, partner_interface::Interface, partner_node::CategoricalConstraint) = true

breakerParameters(interface::Interface, partner_interface::Interface, partner_node::CategoricalConstraint) = (Message{Categorical, Univariate}, size(partner_node.p))

# Average energy functional
averageEnergy(::Type{CategoricalConstraint}, marg_out::ProbabilityDistribution) = differentialEntropy(marg_out)