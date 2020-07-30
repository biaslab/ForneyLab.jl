export ExpectationConstraint

"""
Description:

    Constraints the marginal of the connected variable to an 
    expectation ∫q(x)g(x)dx = G. The parameter η in the node function 
    is actively adapted s.t. the marginal respects the above constraint.
    
    f(out) = exp(η g(out))

Interfaces:

    1. out

Construction:

    ExpectationConstraint(out; g=g, G=G, id=:my_node)
"""
mutable struct ExpectationConstraint <: SoftFactor # TODO: how to handle free energy evaluation for an expectation constrained variable?
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol, Interface}

    g::Function # Constraint function
    G::Float64 # Constraint value
    eta_init::Float64 # Initial value for eta optimization

    function ExpectationConstraint(out; g::Function, G::Float64, eta_init::Float64, id=ForneyLab.generateId(ExpectationConstraint))
        @ensureVariables(out)
        self = new(id, Vector{Interface}(undef, 1), Dict{Symbol,Interface}(), g, G, eta_init)
        ForneyLab.addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)

        return self
    end
end    

slug(::Type{ExpectationConstraint}) = "E"

# A breaker message is required if interface is partnered with an expectation constraint
requiresBreaker(interface::Interface, partner_interface::Interface, partner_node::ExpectationConstraint) = true

breakerParameters(interface::Interface, partner_interface::Interface, partner_node::ExpectationConstraint) = (Message{GaussianMeanVariance, Univariate}, ()) # Univariate only