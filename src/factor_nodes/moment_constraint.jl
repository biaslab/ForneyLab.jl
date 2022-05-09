export MomentConstraint

"""
Description:

    Constraints the marginal of the connected variable to an 
    expectation ∫q(x)g(x)dx = G. The parameter η in the node function 
    is actively adapted s.t. the marginal respects the above constraint.
    Implementation according to (van de Laar et al. "Chance-Constrained
    Active Inference", MIT Neural Computation, 2021).

    f(out) = exp(η g(out))

Interfaces:

    1. out

Construction:

    MomentConstraint(out; g=g, G=G, id=:my_node)
"""
mutable struct MomentConstraint <: SoftFactor
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol, Interface}

    g::Function # Constraint function
    G::Float64 # Constraint value
    eta_init::Float64 # Initial value for eta optimization

    function MomentConstraint(out; g::Function, G::Float64, eta_init::Float64, id=ForneyLab.generateId(MomentConstraint))
        @ensureVariables(out)
        self = new(id, Vector{Interface}(undef, 1), Dict{Symbol,Interface}(), g, G, eta_init)
        ForneyLab.addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)

        return self
    end
end    

slug(::Type{MomentConstraint}) = "E"

# A breaker message is required if interface is partnered with a moment constraint
requiresBreaker(interface::Interface, partner_interface::Interface, partner_node::MomentConstraint) = true

breakerParameters(interface::Interface, partner_interface::Interface, partner_node::MomentConstraint) = (Message{Gaussian{Moments}, Univariate}, ()) # Univariate only

# Constraints do not contribute to average energy
averageEnergy(::Type{MomentConstraint}, marg_out::Distribution) = 0.0
