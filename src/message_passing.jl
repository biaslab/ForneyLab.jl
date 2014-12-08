export  calculateMessage!,
        calculateMessages!,
        calculateForwardMessage!,
        calculateBackwardMessage!,
        executeSchedule,
        clearMessages!

function calculateMessage!(outbound_interface::Interface, graph::FactorGraph=getCurrentGraph())
    # Calculate the outbound message on a specific interface by generating a schedule and executing it.
    # The resulting message is stored in the specified interface and returned.

    # Generate a message passing schedule
    printVerbose("Auto-generating message passing schedule...")
    schedule = generateSchedule!(outbound_interface, graph)
    if verbose show(schedule) end

    # Execute the schedule
    printVerbose("Executing above schedule...")
    executeSchedule(schedule)
    printVerbose("calculateMessage!() done.")

    return outbound_interface.message
end

function pushRequiredInbound!(graph::FactorGraph, inbound_array::Array{Any,1}, node::Node, inbound_interface::Interface, outbound_interface::Interface)
    # Push the inbound message or marginal on inbound_interface, depending on the local graph structure.
    
    if !haskey(graph.edge_to_subgraph, inbound_interface.edge) || !haskey(graph.edge_to_subgraph, outbound_interface.edge)
        # Inbound and/or outbound edge is not explicitly listed in the graph.
        # This is possible if one of those edges is internal to a composite node.
        # We will default to sum-product message passing, and consume the message on the inbound interface.
        # Composite nodes with explicit message passing will throw an error when one of their external interfaces belongs to a different subgraph, so it is safe to assume sum-product.
        try
            return push!(inbound_array, inbound_interface.partner.message)
        catch
            error("$(inbound_interface) is not connected to an edge.")
        end
    end        
    
    inbound_subgraph = graph.edge_to_subgraph[inbound_interface.edge]
    outbound_subgraph = graph.edge_to_subgraph[outbound_interface.edge]

    # Should we require the inbound message or marginal?
    if is(inbound_subgraph, outbound_subgraph)
        # Both edges in same subgraph, require message
        try
            push!(inbound_array, inbound_interface.partner.message)
        catch
            error("$(inbound_interface) is not connected to an edge.")
        end
    else 
        # A subgraph border is crossed, require marginal
        try
            push!(inbound_array, graph.approximate_marginals[(node, inbound_subgraph)])
        catch
            error("Missing approximate marginal for $(inbound_interface)")
        end
    end

end

function updateNodeMessage!(outbound_interface::Interface, graph::FactorGraph=getCurrentGraph())
    # Calculate the outbound message based on the inbound messages and the node update function.
    # The resulting message is stored in the specified interface and is returned.

    node = outbound_interface.node
    # inbound_array holds the inbound messages or marginals on every interface of the node (indexed by the interface id)
    inbound_array = Array(Any, 0)
    outbound_interface_id = 0
    for interface_id = 1:length(node.interfaces)
        interface = node.interfaces[interface_id]
        if interface == outbound_interface
            outbound_interface_id = interface_id
        end
        if (!isdefined(outbound_interface, :dependencies) && outbound_interface_id==interface_id) ||
            (isdefined(outbound_interface, :dependencies) && !(interface in outbound_interface.dependencies))
            # Ignore this interface
            push!(inbound_array, nothing)
        else 
            # Inbound message or marginal is required
            # Put the required inbound message or edge/node marginal for the inbound interface in inbound_array
            pushRequiredInbound!(graph, inbound_array, node, interface, outbound_interface)
        end
    end

    # Evaluate node update function
    printVerbose("Calculate outbound message on $(typeof(node)) $(node.name) interface $outbound_interface_id:")

    return updateNodeMessage!(node, outbound_interface_id, inbound_array...)
end

function calculateMessages!(node::Node)
    # Calculate the outbound messages on all interfaces of node.
    for interface in node.interfaces
        calculateMessage!(interface)
    end
end

# Calculate forward/backward messages on an Edge
calculateForwardMessage!(edge::Edge) = calculateMessage!(edge.tail)
calculateBackwardMessage!(edge::Edge) = calculateMessage!(edge.head)

# Execute schedules
function executeSchedule(schedule::Schedule, graph::FactorGraph=getCurrentGraph())
    # Execute a message passing schedule
    !isempty(schedule) || error("Cannot execute an empty schedule")
    for interface in schedule
        updateNodeMessage!(interface, graph)
    end
    # Return the last message in the schedule
    return schedule[end].message
end
function executeSchedule(schedule::ExternalSchedule, subgraph::Subgraph, graph::FactorGraph=getCurrentGraph())
    # Execute a marginal update schedule
    for entry in schedule
        calculateMarginal!(entry, subgraph, graph)
    end
end
function executeSchedule(subgraph::Subgraph, graph::FactorGraph=getCurrentGraph())
    executeSchedule(subgraph.internal_schedule)
    executeSchedule(subgraph.external_schedule, subgraph, graph)
end

function executeSchedule(graph::FactorGraph=getCurrentGraph())
    for subgraph in graph.factorization
        executeSchedule(subgraph, graph)
    end
end

function clearMessages!(node::Node)
    # Clear all outbound messages on the interfaces of node
    for interface in node.interfaces
        interface.message = nothing
    end
end

function clearMessages!(edge::Edge)
    # Clear all messages on an edge.
    edge.head.message = nothing
    edge.tail.message = nothing
end