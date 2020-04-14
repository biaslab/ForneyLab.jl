export ChanceConstraint

"""
Description:

    Chance constraint on the marginal q of the connected variable x, as
    ϵ ⩽ 1 - ∫_G q(x) dx, where G indicates a region. In other words,
    the probability mass of q is not allowed to overflow the region G
    by more than ϵ.
    
Interfaces:

    1. out

Construction:

    ChanceConstraint(out; G=(min,max), epsilon=epsilon, id=:my_node)
"""
mutable struct ChanceConstraint <: SoftFactor # TODO: how to handle free energy evaluation for a chance-constrained variable?
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol, Interface}

    G::Tuple # Range for G (1-D only)
    epsilon::Float64 # Threshold

    function ChanceConstraint(out; G::Tuple, epsilon::Float64, id=ForneyLab.generateId(ChanceConstraint))
        @ensureVariables(out)
        self = new(id, Vector{Interface}(undef, 1), Dict{Symbol,Interface}(), G, epsilon)
        ForneyLab.addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)

        return self
    end
end    

slug(::Type{ChanceConstraint}) = "P"
