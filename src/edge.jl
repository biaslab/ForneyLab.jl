export Edge
export setForwardMessage!, setBackwardMessage!, getForwardMessage, getBackwardMessage

type Edge <: AbstractEdge
    # An Edge joins two interfaces and has a direction (from tail to head).
    # Edges are mostly useful for code readability, they are not used internally.
    # Forward messages flow in the direction of the Edge (tail to head).
    # Edges can contain marginals, which are the product of the forward and backward message.

    tail::Interface
    head::Interface
    marginal::Any

    function Edge(tail::Interface, head::Interface, forward_message_payload_type::DataType, backward_message_payload_type::DataType; add_to_graph::Bool=true)
        # add_to_graph is false for edges that are internal in a composite node

        (forward_message_payload_type <: MessagePayload) || error("Forward message payload type $(forward_message_payload_type) not supported")
        (backward_message_payload_type <: MessagePayload) || error("Backward message payload type $(backward_message_payload_type) not supported")
        if tail.message != nothing && (typeof(tail.message.payload) != forward_message_payload_type)
            error("Existing forward message ($(typeof(tail.message))) does not match expected type Message{$(forward_message_payload_type)}")
        end
        if head.message != nothing && (typeof(head.message.payload) != backward_message_payload_type)
            error("Existing backward message ($(typeof(head.message))) does not match expected type Message{$(backward_message_payload_type)}")
        end
        (!is(head.node, tail.node)) || error("Cannot connect two interfaces of the same node: ", typeof(head.node), " ", head.node.name)
        # Reset connection between possible previous partners
        if head.partner != nothing
            head.partner.partner = nothing
            head.partner = nothing
        end 
        if tail.partner != nothing
            tail.partner.partner = nothing
            tail.partner = nothing
        end 

        self = new(tail, head, nothing)
        # Assign pointed to edge from interfaces
        tail.edge = self
        head.edge = self
        # Partner head and tail, and merge their families
        tail.partner = head
        head.partner = tail
        # Set expected outbound interface message types
        tail.message_payload_type = forward_message_payload_type
        head.message_payload_type = backward_message_payload_type

        # Backreferences for tail's children
        child_interface = tail.child
        while child_interface != nothing
            child_interface.partner = tail.partner
            child_interface.edge = self
            child_interface.message_payload_type = forward_message_payload_type
            child_interface = child_interface.child
        end
        # Backreferences for head's children
        child_interface = head.child
        while child_interface != nothing
            child_interface.partner = head.partner
            child_interface.edge = self
            child_interface.message_payload_type = backward_message_payload_type
            child_interface = child_interface.child
        end

        # Incorporate edge and nodes in current graph
        if add_to_graph
            graph = getCurrentGraph()
            (length(graph.factorization) == 1) || error("Cannot create Edge in an already factorized graph; first build the graph, then define factorizations.")
            subgraph = graph.factorization[1] # There is only one
            graph.edge_to_subgraph[self] = subgraph # Add edge to internal mapping
            push!(subgraph.internal_edges, self) # Define edge as internal
            push!(subgraph.nodes, tail.node) # Add node to subgraph
            push!(subgraph.nodes, head.node)
        end

        return self
    end
end
Edge(tail::Interface, head::Interface, message_payload_type::DataType; args...) = Edge(tail, head, message_payload_type, message_payload_type; args...)
function Edge(tail::Interface, head::Interface; args...)
    forward_message_payload_type = backward_message_payload_type = GaussianDistribution # Default to Gaussian; we could also do a check whether the nodes accept Gaussians
    (tail.message == nothing) || (forward_message_payload_type = typeof(tail.message.payload))
    (head.message == nothing) || (backward_message_payload_type = typeof(head.message.payload))
    Edge(tail, head, forward_message_payload_type, backward_message_payload_type; args...)
end

# Edge constructors that accept nodes instead of a specific Interface
# firstFreeInterface(node) should be overloaded for nodes with interface-invariant node functions
firstFreeInterface(node::Node) = error("Cannot automatically pick a free interface on non-symmetrical $(typeof(node)) $(node.name)")
Edge(tail_node::Node, head::Interface; args...) = Edge(firstFreeInterface(tail_node), head; args...)
Edge(tail_node::Node, head::Interface, message_payload_type::DataType; args...) = Edge(firstFreeInterface(tail_node), head, message_payload_type; args...)
Edge(tail_node::Node, head::Interface, forward_message_payload_type::DataType, backward_message_payload_type::DataType; args...) = Edge(firstFreeInterface(tail_node), head, forward_message_payload_type, backward_message_payload_type; args...)

Edge(tail::Interface, head_node::Node; args...) = Edge(tail, firstFreeInterface(head_node); args...)
Edge(tail::Interface, head_node::Node, message_payload_type::DataType; args...) = Edge(tail, firstFreeInterface(head_node), message_payload_type; args...)
Edge(tail::Interface, head_node::Node, forward_message_payload_type::DataType, backward_message_payload_type::DataType; args...) = Edge(tail, firstFreeInterface(head_node), forward_message_payload_type, backward_message_payload_type; args...)

Edge(tail_node::Node, head_node::Node; args...) = Edge(firstFreeInterface(tail_node), firstFreeInterface(head_node); args...)
Edge(tail_node::Node, head_node::Node, message_payload_type::DataType; args...) = Edge(firstFreeInterface(tail_node), firstFreeInterface(head_node), message_payload_type; args...)
Edge(tail_node::Node, head_node::Node, forward_message_payload_type::DataType, backward_message_payload_type::DataType; args...) = Edge(firstFreeInterface(tail_node), firstFreeInterface(head_node), forward_message_payload_type, backward_message_payload_type; args...)

function show(io::IO, edge::Edge)
    println(io, "Edge from $(typeof(edge.tail.node)) $(edge.tail.node.name):$(findfirst(edge.tail.node.interfaces, edge.tail)) to $(typeof(edge.head.node)) $(edge.head.node.name):$(findfirst(edge.head.node.interfaces, edge.head)).")
    println(io, "Accepted forward message payload type: $(edge.tail.message_payload_type).")
    println(io, "Accepted backward message payload type: $(edge.head.message_payload_type).")
end

setForwardMessage!(edge::Edge, message::Message) = setMessage!(edge.tail, message)
setBackwardMessage!(edge::Edge, message::Message) = setMessage!(edge.head, message)
getForwardMessage(edge::Edge) = edge.tail.message
getBackwardMessage(edge::Edge) = edge.head.message

function getOrCreateMarginal(edge::Edge, assign_distribution::DataType)
    # Looks for a marginal on edge.
    # When no marginal is present, it sets and returns an uninformative distribution.
    # Otherwise it returns the present marginal. Used for fast marginal calculations.
    if edge.marginal==nothing
        if assign_distribution <: ProbabilityDistribution 
            edge.marginal = uninformative(assign_distribution)
        else
            error("Unknown assign type argument")
        end
    end
    return edge.marginal
end
