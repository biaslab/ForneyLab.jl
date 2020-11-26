export BernoulliConstraint

"""
Description:

    Constraints the marginal of the connected variable to a Bernoulli belief.
    
Interfaces:

    1. out

Construction:

    BernoulliConstraint(out; p=p, id=:my_node)
"""
mutable struct BernoulliConstraint <: SoftFactor # TODO: how to handle free energy evaluation for a point-mass constrained variable?
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol, Interface}

    p::Float64

    function BernoulliConstraint(out; p::Float64,  id=ForneyLab.generateId(BernoulliConstraint))
        @ensureVariables(out)
        self = new(id, Vector{Interface}(undef, 1), Dict{Symbol,Interface}(), p)
        ForneyLab.addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)

        return self
    end
end    

slug(::Type{BernoulliConstraint}) = "Ber_q"

# A breaker message is required if interface is partnered with a point-mass constraint
requiresBreaker(interface::Interface, partner_interface::Interface, partner_node::BernoulliConstraint) = true

breakerParameters(interface::Interface, partner_interface::Interface, partner_node::BernoulliConstraint) = (Message{Bernoulli, Univariate}, ())

# Average energy functional
averageEnergy(::Type{BernoulliConstraint}, marg_out::ProbabilityDistribution) = differentialEntropy(marg_out)