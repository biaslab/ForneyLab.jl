export Edge

"""
An Edge joins two interfaces (half-edges) `a` and `b`.
"""
type Edge <: AbstractEdge
    id::Symbol
    a::Interface
    b::Interface

    function Edge(a::Interface, b::Interface; id=Symbol("$(a.node.id)_$(b.node.id)"))
        # add_to_graph is false for edges that are internal in a composite node
        current_graph = currentGraph()
        (b.partner == nothing && a.partner == nothing) || error("Previously defined edges cannot be repositioned.")
        (hasNode(current_graph, b.node) && hasNode(current_graph, a.node)) || error("Node does not belong to the current graph.")
        !haskey(current_graph.edges, id) || error("The edge id $(id) already exists in the current graph. Consider specifying an explicit id.")

        self = new(id, a, b)

        # Assign pointed to edge from interfaces
        a.edge = self
        b.edge = self
        # Partner b and a, and merge their families
        a.partner = b
        b.partner = a

        # Add edge to current_graph
        current_graph.edges[self.id] = self
        current_graph.prepared_algorithm = nothing # Modifying the graph 'unprepares' any InferenceAlgorithm

        return self
    end
end

function show(io::IO, edge::Edge)
    if (a_handle = handle(edge.a)) != ""
        a_interface = ((typeof(a_handle)==Symbol) ? "i[:$(a_handle)]" : "i[$(a_handle)]")
    else
        a_interface = "interfaces[$(findfirst(edge.a.node.interfaces, edge.a))]"
    end
    if (b_handle = handle(edge.b)) != ""
        b_interface = ((typeof(b_handle)==Symbol) ? "i[:$(b_handle)]" : "i[$(b_handle)]")
    else
        b_interface = "interfaces[$(findfirst(edge.b.node.interfaces, edge.b))]"
    end
    println(io, "Edge with id $(edge.id) from $(edge.a.node.id).$(a_interface) to $(edge.b.node.id).$(b_interface).")
end

function show(io::IO, edges::Union{Vector{Edge}, Set{Edge}})
    println(io, "Edges:")
    for edge in edges
        show(io, edge)
    end
end

"""Compare edges by alphabetically comparing the ids in the order (b, a)"""
Base.isless(e1::Edge, e2::Edge) = isless("$(e1.id)", "$(e2.id)")
Base.isless(e_arr1::Array{Edge}, e_arr2::Array{Edge}) = isless("$(e_arr1[1,1].id)", "$(e_arr2[1,1].id)")
Base.isless(e1::Edge, e_arr2::Array{Edge}) = isless("$(e1.id)", "$(e_arr2[1,1].id)")
Base.isless(e_arr1::Array{Edge}, e2::Edge) = isless("$(e_arr1[1,1].id)", "$(e2.id)")