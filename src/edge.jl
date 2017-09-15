export Edge

"""
An Edge joins two interfaces (half-edges) `a` and `b`.
"""
type Edge <: AbstractEdge
    variable::AbstractVariable
    a::Union{Interface, Void}
    b::Union{Interface, Void}

    function Edge(var::AbstractVariable, a::Interface)
        current_graph = currentGraph()
        (a.partner == nothing && a.edge == nothing) || error("Interface is already connected to an Edge.")
        hasNode(current_graph, a.node) || error("Node does not belong to the current graph.")
        hasVariable(current_graph, var) || error("Variable does not belong to the current graph.")

        self = new(var, a, nothing)
        a.edge = self

        push!(var.edges, self)
        push!(current_graph.edges, self)

        return self
    end
end

Edge(var::AbstractVariable, a::Interface, b::Interface) = connect!(Edge(var, a), b)

"""
Connect loose end of edge to interface b.
"""
function connect!(edge::Edge, b::Interface)
    (edge.b == nothing) || error("The edge does not have a loose end.")
    (b.partner == nothing && b.edge == nothing) || error("Interface is already connected to an Edge.")
    hasNode(currentGraph(), b.node) || error("Node does not belong to the current graph.")
    
    edge.b = b
    b.edge = edge
    edge.a.partner = b
    b.partner = edge.a

    return edge
end

"""
Disconnect edge from interface.
"""
function disconnect!(edge::Edge, interface::Interface)
    if !is(edge.a, interface) && !is(edge.b, interface)
        error("Cannot disconnect from an interface that is not connected.")
    end

    edge.a.partner = nothing
    edge.b.partner = nothing
    if is(edge.a, interface)
        edge.a.edge = nothing
        edge.a = edge.b
    else 
        edge.b.edge = nothing
    end
    edge.b = nothing

    return edge
end

function show(io::IO, edge::Edge)
    if edge.a != nothing
        a = "$(edge.a.node.id)."
        if (a_handle = handle(edge.a)) != ""
            a *= ((typeof(a_handle)==Symbol) ? "i[:$(a_handle)]" : "i[$(a_handle)]")
        else
            a *= "interfaces[$(findfirst(edge.a.node.interfaces, edge.a))]"
        end
    else
        a = "NONE"
    end

    if edge.b != nothing
        b = "$(edge.b.node.id)."
        if (b_handle = handle(edge.b)) != ""
            b *= ((typeof(b_handle)==Symbol) ? "i[:$(b_handle)]" : "i[$(b_handle)]")
        else
            b *= "interfaces[$(findfirst(edge.b.node.interfaces, edge.b))]"
        end
    else
        b = "NONE"
    end   

    println(io, "Edge belonging to variable $(edge.variable.id): ( $(a) )----( $(b) ).")
end

Base.isless(e1::Edge, e2::Edge) = isless(e1.variable, e2.variable)

function show(io::IO, edges::Union{Vector{Edge}, Set{Edge}})
    println(io, "Edges:")
    for edge in edges
        show(io, edge)
    end
end
