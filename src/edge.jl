export Edge
export setForwardMessage!, setBackwardMessage!, forwardMessage, backwardMessage, ensureMarginal!

"""
An Edge joins two interfaces and has a direction (from tail to head).
Forward messages flow in the direction of the Edge (tail to head).
An Edge can contain a marginal, which is the product of the forward and backward message.
"""
type Edge <: AbstractEdge
    id::Symbol
    tail::Interface
    head::Interface
    marginal::Union{ProbabilityDistribution, Void}

    function Edge(tail::Interface, head::Interface; id=symbol("$(tail.node.id)_$(head.node.id)"))
        # add_to_graph is false for edges that are internal in a composite node
        current_graph = currentGraph()
        !is(head.node, tail.node) || error("Cannot connect two interfaces of the same node: $(typeof(head.node)) $(head.node.id)")
        (head.partner == nothing && tail.partner == nothing) || error("Previously defined edges cannot be repositioned.")
        !current_graph.locked || error("Cannot extend a locked FactorGraph.")
        hasNode(current_graph, head.node) || error("Head node does not belong to the current graph.")
        hasNode(current_graph, tail.node) || error("Tail node does not belong to the current graph.")
        !haskey(current_graph.edges, id) || error("The edge id $(id) already exists in the current graph. Consider specifying an explicit id.")

        self = new(id, tail, head, nothing)

        # Assign pointed to edge from interfaces
        tail.edge = self
        head.edge = self
        # Partner head and tail, and merge their families
        tail.partner = head
        head.partner = tail

        # Add edge to current_graph
        current_graph.edges[self.id] = self

        return self
    end
end

Base.deepcopy(::Edge) = error("deepcopy(::Edge) is not possible. You should construct a new Edge, or copy the entire FactorGraph using deepcopy(::FactorGraph).")

# Edge constructors that accept nodes instead of a specific Interface
# firstFreeInterface(node) should be overloaded for nodes with interface-invariant node functions
firstFreeInterface(node::Node) = error("Cannot automatically pick a free interface on non-symmetrical $(typeof(node)) $(node.id)")
Edge(tail_node::Node, head::Interface; args...) = Edge(firstFreeInterface(tail_node), head; args...)
Edge(tail::Interface, head_node::Node; args...) = Edge(tail, firstFreeInterface(head_node); args...)
Edge(tail_node::Node, head_node::Node; args...) = Edge(firstFreeInterface(tail_node), firstFreeInterface(head_node); args...)


function show(io::IO, edge::Edge)
    if (tail_handle = handle(edge.tail)) != ""
        tail_interface = ((typeof(tail_handle)==Symbol) ? "i[:$(tail_handle)]" : "i[$(tail_handle)]")
    else
        tail_interface = "interfaces[$(findfirst(edge.tail.node.interfaces, edge.tail))]"
    end
    if (head_handle = handle(edge.head)) != ""
        head_interface = ((typeof(head_handle)==Symbol) ? "i[:$(head_handle)]" : "i[$(head_handle)]")
    else
        head_interface = "interfaces[$(findfirst(edge.head.node.interfaces, edge.head))]"
    end
    println(io, "Edge with id $(edge.id) from $(edge.tail.node.id).$(tail_interface) to $(edge.head.node.id).$(head_interface).")
end

function show(io::IO, edges::Union{Vector{Edge}, Set{Edge}})
    println(io, "Edges:")
    for edge in edges
        show(io, edge)
    end
end

setForwardMessage!(edge::Edge, message::Message) = setMessage!(edge.tail, message)
setBackwardMessage!(edge::Edge, message::Message) = setMessage!(edge.head, message)
forwardMessage(edge::Edge) = edge.tail.message
backwardMessage(edge::Edge) = edge.head.message

function ensureMarginal!{T<:ProbabilityDistribution}(edge::Edge, distribution_type::Type{T})
    # Ensure that edge carries a marginal of type distribution_type, used for in place updates
    if edge.marginal==nothing || (typeof(edge.marginal) <: distribution_type)==false
        if distribution_type <: DeltaDistribution{Float64}
            edge.marginal = DeltaDistribution() # vague() not implemented
        else
            edge.marginal = vague(distribution_type)
        end
    end

    return edge.marginal
end

# Compare edges by alphabetically comparing the ids in the order (head, tail)
Base.isless(e1::Edge, e2::Edge) = isless("$(e1.id)", "$(e2.id)")
Base.isless(e_arr1::Array{Edge}, e_arr2::Array{Edge}) = isless("$(e_arr1[1,1].id)", "$(e_arr2[1,1].id)")
Base.isless(e1::Edge, e_arr2::Array{Edge}) = isless("$(e1.id)", "$(e_arr2[1,1].id)")
Base.isless(e_arr1::Array{Edge}, e2::Edge) = isless("$(e_arr1[1,1].id)", "$(e2.id)")
