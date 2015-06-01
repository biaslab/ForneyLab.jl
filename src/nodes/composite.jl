export CompositeNode, addRule!, nodes, isDeterministic

type CompositeNode <: Node
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol,Interface}
    internal_graph::FactorGraph
    computation_rules::Dict{(Interface,Function),Algorithm}
    interfaceid_to_terminalnode::Array{TerminalNode,1}
    terminalnode_to_interface::Dict{TerminalNode,Interface}
    deterministic::Bool
end

function CompositeNode(graph::FactorGraph, terminals...; id=generateNodeId(), deterministic=false)
    # Convert graph into a CompositeNode.
    # terminals... is an array of TerminalNodes in graph, that will be bound to interfaces of the CompositeNode.
    # The ids of the terminals will be used as interface handles.
    # This function creates a new current FactorGraph so the user can continue working in the higher level graph.
    self = CompositeNode(id, Interface[], Dict{Symbol,Interface}(), graph, Dict{(Interface,Function),Algorithm}(), TerminalNode[], Dict{TerminalNode,Interface}(), deterministic)

    for terminal in terminals
        (typeof(terminal) == TerminalNode) || error("Only a TerminalNode can be bound to an Interface of the CompositeNode, not $(typeof(terminal)).")
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
    new_g.n[id] = self # Add newly created composite node to current graph

    return self    
end

isDeterministic(composite_node::CompositeNode) = composite_node.deterministic

function addRule!(composite_node::CompositeNode, outbound_interface::Interface, message_calculation_rule::Function, algorithm::Algorithm)
    (outbound_interface.node == composite_node) || error("The outbound interface does not belong to the specified composite node $(composite_node.id)")
    composite_node.computation_rules[(outbound_interface,message_calculation_rule)] = algorithm

    return composite_node
end

nodes(node::CompositeNode) = nodes(node.internal_graph)

function sumProduct!(node::CompositeNode,
                     outbound_interface_id::Int,
                     inbounds...)
    outbound_interface = node.interfaces[outbound_interface_id]
    internal_outbound_interface = node.interfaceid_to_terminalnode[outbound_interface_id].interfaces[1].partner

    # Check if there is a computation rule defined for this case
    if !haskey(node.computation_rules, (outbound_interface, sumProduct!))
        # Try to automatically generate a sum-product algorithm
        clearMessages!(node.internal_graph)
        algo = SumProduct.Algorithm(internal_outbound_interface)
        node.computation_rules[(outbound_interface, sumProduct!)] = algo
    end

    # Move all inbound messages to corresponding terminal nodes in the internal graph
    for interface_id=1:length(node.interfaces)
        if interface_id != outbound_interface_id 
            node.interfaceid_to_terminalnode[interface_id].value = node.interfaces[interface_id].partner.message.payload
        end
    end

    # Execute the correct algorithm
    run(node.computation_rules[(outbound_interface, sumProduct!)], node.internal_graph)

    # Move the calculated messages to the interfaces of node
    outbound_interface.message = internal_outbound_interface.message

    return(:empty, outbound_interface.message)
end
