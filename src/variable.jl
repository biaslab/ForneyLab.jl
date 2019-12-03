export
Variable

"""
A `Variable` encompasses one or more edges in a `FactorGraph`.
"""
mutable struct Variable <: AbstractVariable
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
        Edge(var, iface)
    elseif var.edges[end].b == nothing
        # Connect to the loose end of an existing Edge
        connect!(var.edges[end], iface)
    else
        # Insert an equality constraint node
        equ_idx = Int((length(var.edges) - 1) / 2) + 1
        equ = Equality(id=Symbol("equ_$(var.id)_$(equ_idx)"))
        disconnected_iface = var.edges[end].b
        disconnect!(var.edges[end], disconnected_iface)
        connect!(var.edges[end], equ.interfaces[1])
        Edge(var, disconnected_iface, equ.interfaces[2])
        Edge(var, iface, equ.interfaces[3])
    end

    return iface
end

"""
Collect all edges corresponding with variable(s)
"""
edges(variable::Variable) = Set{Edge}(variable.edges)
edges(variables::Set{Variable}) = union(Set((edges(v) for v=variables))...)

Base.isless(v1::Variable, v2::Variable) = isless("$(v1.id)", "$(v2.id)")