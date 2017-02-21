export Variable

"""
A Variable encompasses one or more edges in a FactorGraph.
"""
immutable Variable <: AbstractVariable
    id::Symbol
    edges::Vector{Edge}

    function Variable(;id=generateId(Variable))
        self = new(id, Vector{Edge}())
        addVariable!(currentGraph(), self)
        return self
    end
end

"""
`associate!(interface, variable)` associates `interface` with
`variable` by connecting `interface` to an `Edge` belonging to `variable`.
"""
function associate!(iface::Interface, var::Variable)
    if isempty(var.edges)
        # Make a new Edge
        edge = Edge(iface)
    elseif var.edges[end].b == nothing
        # Connect to the loose end of an existing Edge

        # update partner, update edge, update interface
    else
        # Insert an equality constraint node

    end

    return iface
end