export CompositeNode, addRule!, nodes, isDeterministic

type CompositeNode <: Node
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol,Interface}
    internal_graph::FactorGraph
    computation_rules::Dict{Interface, InferenceAlgorithm} # Maps external interfaces to inference algorithm
    interfaceid_to_terminalnode::Array{TerminalNode,1} # Maps external interface ids to internal terminal nodes
    terminalnode_to_interface::Dict{TerminalNode,Interface} # Maps internal terminal nodes to external interfaces
    deterministic::Bool
end

function CompositeNode(graph::FactorGraph, terminals...; id=generateNodeId(CompositeNode), deterministic=false)
    # Convert graph into a CompositeNode.
    # terminals... is an array of TerminalNodes in graph, that will be bound to interfaces of the CompositeNode.
    # The ids of the terminals will be used as interface handles.
    # This function creates a new current FactorGraph so the user can continue working in the higher level graph.
    self = CompositeNode(id, Interface[], Dict{Symbol,Interface}(), graph, Dict{Interface, InferenceAlgorithm}(), TerminalNode[], Dict{TerminalNode,Interface}(), deterministic)

    for terminal in terminals
        hasNode(graph, terminal) || error("$(node.id) not in graph")
        (typeof(terminal) <: TerminalNode) || error("Only a TerminalNode can be bound to an Interface of the CompositeNode, not $(typeof(terminal)).")
        if terminal in self.interfaceid_to_terminalnode
            error("Cannot bind the same TerminalNode to multiple interfaces")
        end
        push!(self.interfaces, Interface(self))
        self.i[terminal.id] = self.interfaces[end]
        push!(self.interfaceid_to_terminalnode, terminal)
        self.terminalnode_to_interface[terminal] = self.interfaces[end]
    end
    self.internal_graph.locked = true
    new_g = FactorGraph()
    addNode!(new_g, self) # Add newly created composite node to new current graph

    return self
end
CompositeNode(terminals...; id=generateNodeId(CompositeNode), deterministic=false) = CompositeNode(current_graph, terminals..., id=id, deterministic=deterministic)

isDeterministic(composite_node::CompositeNode) = composite_node.deterministic

function addRule!(composite_node::CompositeNode, outbound_interface::Interface, algorithm::InferenceAlgorithm)
    (outbound_interface.node == composite_node) || error("The outbound interface does not belong to the specified composite node $(composite_node.id)")
    composite_node.computation_rules[outbound_interface] = algorithm

    return composite_node
end

nodes(node::CompositeNode) = nodes(node.internal_graph)

function sumProduct!{   T<:Val}(node::CompositeNode,
                        outbound_interface_id::Type{T},
                        args...)

    outbound_dist = args[end] # Dummy declaration, just to be clear on which one the outbound distribution is 
    outbound_interface_id = outbound_interface_id.parameters[1]
    outbound_interface = node.interfaces[outbound_interface_id]

    # Move all inbound messages to corresponding terminal nodes in the internal graph
    for (id, interface) in enumerate(node.interfaces)
        if id != outbound_interface_id
            node.interfaceid_to_terminalnode[id].value = interface.partner.message.payload
        end
    end

    # Execute the internal graph algorithm
    parent_graph = currentGraph()
    setCurrentGraph(node.internal_graph)
    run(node.computation_rules[outbound_interface]) # Executes an in-place operation on outbound_dist
    # reset the graph
    setCurrentGraph(parent_graph)

    # Move the calculated messages from the internal to the external interface
    internal_outbound_interface = node.interfaceid_to_terminalnode[outbound_interface_id].interfaces[1].partner
    outbound_interface.message = internal_outbound_interface.message

    return outbound_dist # Return the altered outbound_dist
end
