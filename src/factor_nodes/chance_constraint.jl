export ChanceConstraint

"""
Description:

    Chance constraint on the marginal q of the connected variable x, as
    ϵ ⩽ 1 - ∫_G q(x) dx, where g(x) ∈ {0, 1} indicates whether x ∈ G.
    
Interfaces:

    1. out

Construction:

    ChanceConstraint(out; g=g, epsilon=epsilon, id=:my_node)
"""
mutable struct ChanceConstraint <: SoftFactor # TODO: how to handle free energy evaluation for a chance-constrained variable?
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol, Interface}

    g::Function # Indicator function for x ∈ G
    epsilon::Float64 # Threshold

    function ChanceConstraint(out; g::Function, epsilon::Float64, id=ForneyLab.generateId(ChanceConstraint))
        @ensureVariables(out)
        self = new(id, Vector{Interface}(undef, 1), Dict{Symbol,Interface}(), g, epsilon)
        ForneyLab.addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)

        return self
    end
end    

slug(::Type{ChanceConstraint}) = "P"
