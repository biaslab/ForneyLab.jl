export  calculateMessage!,
        calculateForwardMessage!,
        calculateBackwardMessage!,
        execute,
        clearMessages!

function calculateMessage!(outbound_interface::Interface, graph::FactorGraph=currentGraph())
    # Calculate the outbound message on a specific interface by generating a schedule and executing it.
    # The resulting message is stored in the specified interface and returned.

    # Generate a message passing schedule
    printVerbose("Auto-generating message passing schedule...")
    schedule = generateSchedule!(outbound_interface, graph)
    if verbose show(schedule) end

    # Execute the schedule
    printVerbose("Executing above schedule...")
    execute(schedule)
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

function execute(schedule_entry::ScheduleEntry, graph::FactorGraph=currentGraph())
    # Calculate the outbound message based on the inbound messages and the message calculation rule.
    # The resulting message is stored in the specified interface and is returned.

    outbound_interface = schedule_entry.interface

    # Preprocessing: collect all inbound messages
    node = outbound_interface.node
    inbound_array = Array(Any, 0) # inbound_array holds the inbound messages or marginals on every interface of the node (indexed by the interface id)
    outbound_interface_id = 0
    for interface_id = 1:length(node.interfaces)
        interface = node.interfaces[interface_id]
        if interface == outbound_interface
            outbound_interface_id = interface_id
        end
        if (outbound_interface_id==interface_id)
            # Ignore this interface
            # In the future we might want to have decent dependency checking here
            push!(inbound_array, nothing)
        else
            # Inbound message or marginal is required
            # Put the required inbound message or edge/node marginal for the inbound interface in inbound_array
            pushRequiredInbound!(graph, inbound_array, node, interface, outbound_interface)
        end
    end

    # Evaluate message calculation rule
    printVerbose("Calculate outbound message on $(typeof(node)) $(node.name) interface $outbound_interface_id:")
    outbound_message = schedule_entry.message_calculation_rule(node, outbound_interface_id, inbound_array...)

    # Post processing?
    if isdefined(schedule_entry, :post_processing)
        outbound_message = node.interfaces[outbound_interface_id].message = Message(schedule_entry.post_processing(outbound_message.payload))
    end

    return outbound_message
end

# Calculate forward/backward messages on an Edge
calculateForwardMessage!(edge::Edge) = calculateMessage!(edge.tail)
calculateBackwardMessage!(edge::Edge) = calculateMessage!(edge.head)

# Execute schedules
function execute(schedule::Any, graph::FactorGraph=currentGraph())
    # Execute a message passing schedule
    !isempty(schedule) || error("Cannot execute an empty schedule")
    for schedule_entry in schedule
        execute(schedule_entry, graph)
    end
    # Return the last message in the schedule
    return schedule[end].interface.message
end
function execute(schedule::ExternalSchedule, subgraph::Subgraph, graph::FactorGraph=currentGraph())
    # Execute a marginal update schedule
    for entry in schedule
        calculateMarginal!(entry, subgraph, graph)
    end
end
function execute(subgraph::Subgraph, graph::FactorGraph=currentGraph())
    execute(subgraph.internal_schedule)
    execute(subgraph.external_schedule, subgraph, graph)
end

function execute(graph::FactorGraph=currentGraph())
    for subgraph in graph.factorization
        execute(subgraph, graph)
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
