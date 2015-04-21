function pushRequiredInbound!(scheme::InferenceScheme, inbound_array::Array{Any,1}, node::Node, inbound_interface::Interface, outbound_interface::Interface)
    # Push the inbound message or marginal on inbound_interface, depending on the local graph structure.

    if !haskey(scheme.edge_to_subgraph, inbound_interface.edge) || !haskey(scheme.edge_to_subgraph, outbound_interface.edge)
        # Inbound and/or outbound edge is not explicitly listed in the scheme.
        # This is possible if one of those edges is internal to a composite node.
        # We will default to sum-product message passing, and consume the message on the inbound interface.
        # Composite nodes with explicit message passing will throw an error when one of their external interfaces belongs to a different subgraph, so it is safe to assume sum-product.
        try return push!(inbound_array, inbound_interface.partner.message) catch error("$(inbound_interface) is not connected to an edge.") end
    end

    # Should we require the inbound message or marginal?
    if is(subgraph(scheme, inbound_interface.edge), subgraph(scheme, outbound_interface.edge))
        # Both edges in same subgraph, require message
        try push!(inbound_array, inbound_interface.partner.message) catch error("$(inbound_interface) is not connected to an edge.") end
    else
        # A subgraph border is crossed, require marginal
        # The factor is the set of internal edges that are in the same subgraph
        try push!(inbound_array, qDistribution(scheme, node, inbound_interface.edge)) catch error("Missing approximate marginal for $(inbound_interface)") end
    end

end

function execute(subgraph::Subgraph)
    printVerbose("Subgraph $(findfirst(scheme.factorization, subgraph)):")
    # Execute internal schedule
    execute(subgraph.internal_schedule)
    # Update q-distributions at external edges
    g_nodes = nodesConnectedToExternalEdges(subgraph)
    for node in g_nodes
        calculateQDistribution!(node, qFactor(scheme, node, subgraph), scheme)
    end
end

function execute(scheme::InferenceScheme)
    for subgraph in scheme.factorization
        execute(subgraph, scheme)
    end
end
