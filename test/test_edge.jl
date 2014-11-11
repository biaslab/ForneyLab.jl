#####################
# Integration tests
#####################

facts("Edge integration tests") do
    context("Edge constructor should add edge to edge_to_subgraph mapping") do
        (node1, node2) = initializePairOfNodes()
        # Couple the interfaces that carry GeneralMessage
        edge = Edge(node2.interfaces[1], node1.interfaces[1]) # Edge from node 2 to node 1
        graph = getCurrentGraph()
        @fact graph.edge_to_subgraph[edge] => graph.factorization[1]
    end

    context("Nodes can be coupled by edges using the explicit interface names") do
        (node1, node2) = initializePairOfNodes()
        # Couple the interfaces that carry GeneralMessage
        edge = Edge(node2.out, node1.in1) # Edge from node 2 to node 1
        testInterfaceConnections(node1, node2)
    end

    context("Edge should throw an error when two interfaces on the same node are connected") do
        node = FixedGainNode()
        # Connect output directly to input
        @fact_throws Edge(node.interfaces[2], node.interfaces[1])
    end

    context("Edge constructor should write the expected message value types to the interfaces") do
        (node1, node2) = initializePairOfMockNodes()
        edge = Edge(node1.out, node2.out, GaussianDistribution, Float64)
        @fact edge.tail.message_payload_type => GaussianDistribution
        @fact edge.head.message_payload_type => Float64
    end

    context("Edge construction should couple interfaces to edge") do
        (node1, node2) = initializePairOfMockNodes()
        @fact node1.out.edge => nothing
        @fact node2.out.edge => nothing
        edge = Edge(node1.out, node2.out)
        @fact node1.out.edge => edge
        @fact node2.out.edge => edge
    end

    context("Edge should couple standard to composite nodes") do
        comp_node = GainEqualityCompositeNode()
        node = TerminalNode()
        edge = Edge(node.out, comp_node.in1)
        @fact comp_node.equality_node.interfaces[1].partner => node.out
        @fact comp_node.equality_node.interfaces[1].edge => edge
    end

    context("Edge should throw an error when the user attempts to reposition") do
        node1 = TerminalNode(name="node1")
        node2 = TerminalNode(name="node2")
        node3 = TerminalNode(name="node3")
        Edge(node1.out, node2.out)
        @fact_throws Edge(node1.out, node3.out)
    end

    context("Edges can be (pseudo) ordered") do
        edge1 = Edge(TerminalNode().out, TerminalNode().out)
        edge2 = Edge(TerminalNode().out, TerminalNode().out)
        @fact edge1 < edge2 => true
    end
end