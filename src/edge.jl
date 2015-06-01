export Edge
export setForwardMessage!, setBackwardMessage!, forwardMessage, backwardMessage, ensureMarginal!

type Edge <: AbstractEdge
    # An Edge joins two interfaces and has a direction (from tail to head).
    # Forward messages flow in the direction of the Edge (tail to head).
    # An Edge can contain a marginal, which is the product of the forward and backward message.
    # If distribution_type is defined, it restricts the distribution type of the marginal.

    id::Symbol
    tail::Interface
    head::Interface
    marginal::Union(ProbabilityDistribution, Nothing)
    distribution_type::DataType

    function Edge(tail::Interface, head::Interface, distribution_type=Any; id=symbol("$(tail.node.id)_$(head.node.id)"))
        # add_to_graph is false for edges that are internal in a composite node
        # Cautionary note: replacing "distribution_type=Any" by "distribution_type::DataType=Any" causes segfaults (?!?)
        (!is(head.node, tail.node)) || error("Cannot connect two interfaces of the same node: $(typeof(head.node)) $(head.node.id)")
        (head.partner == nothing && tail.partner == nothing) || error("Previously defined edges cannot be repositioned.")
        (current_graph.locked == false) || error("Cannot extend a locked FactorGraph.")

        self = new(id, tail, head, nothing, distribution_type)

        # Assign pointed to edge from interfaces
        tail.edge = self
        head.edge = self
        # Partner head and tail, and merge their families
        tail.partner = head
        head.partner = tail

        # Add edge to current_graph
        current_graph.e[self.id] = self

        return self
    end
end

# Edge constructors that accept nodes instead of a specific Interface
# firstFreeInterface(node) should be overloaded for nodes with interface-invariant node functions
firstFreeInterface(node::Node) = error("Cannot automatically pick a free interface on non-symmetrical $(typeof(node)) $(node.id)")
Edge(tail_node::Node, head::Interface, distribution_type=Any; args...) = Edge(firstFreeInterface(tail_node), head, distribution_type; args...)
Edge(tail::Interface, head_node::Node, distribution_type=Any; args...) = Edge(tail, firstFreeInterface(head_node), distribution_type; args...)
Edge(tail_node::Node, head_node::Node, distribution_type=Any; args...) = Edge(firstFreeInterface(tail_node), firstFreeInterface(head_node), distribution_type; args...)


function show(io::IO, edge::Edge)
    println(io, "Edge from $(typeof(edge.tail.node)) $(edge.tail.node.id):$(findfirst(edge.tail.node.interfaces, edge.tail)) to $(typeof(edge.head.node)) $(edge.head.node.id):$(findfirst(edge.head.node.interfaces, edge.head)).")
    (edge.distribution_type == Any) || println(io, "Marginal distribution type: $(edge.distribution_type).")
end

setForwardMessage!(edge::Edge, message::Message) = setMessage!(edge.tail, message)
setBackwardMessage!(edge::Edge, message::Message) = setMessage!(edge.head, message)
forwardMessage(edge::Edge) = edge.tail.message
backwardMessage(edge::Edge) = edge.head.message

function ensureMarginal!{T<:ProbabilityDistribution}(edge::Edge, distribution_type::Type{T}=edge.distribution_type)
    # Ensure that edge carries a marginal of type distribution_type, used for in place updates
    if edge.marginal==nothing || (typeof(edge.marginal) <: distribution_type)==false
        (distribution_type <: edge.distribution_type) || error("Cannot create marginal of type $(distribution_type) since the edge requires a different marginal distribution type. Edge:\n$(edge)")
        try
            # Initialize the marginal as a vague "uninformative" distribution of the specified type
            edge.marginal = vague(distribution_type)
        catch
            # vague() not implemented
            # Instead, initialize the marginal with the default distribution of the specified type
            edge.marginal = distribution_type()
        end
    end

    return edge.marginal
end

# Compare edges by alphabetically comparing the handles in the order (head, tail)
Base.isless(e1::Edge, e2::Edge) = isless("$(e1.id)", "$(e2.id)")