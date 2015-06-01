#####################
# Integration tests
#####################

facts("Edge integration tests") do
    context("Nodes can be coupled by edges using the explicit interface names") do
        (node1, node2) = initializePairOfNodes()
        # Couple the interfaces that carry GeneralMessage
        edge = Edge(node2.i[:out], node1.i[:in]) # Edge from node 2 to node 1
        testInterfaceConnections(node1, node2)
    end

    context("Edge should throw an error when two interfaces on the same node are connected") do
        node = FixedGainNode()
        # Connect output directly to input
        @fact_throws Edge(node.interfaces[2], node.interfaces[1])
    end

    context("Edge constructor should write the distribution type") do
        (node1, node2) = initializePairOfMockNodes()
        edge = Edge(node1.i[:out], node2.i[:out], GaussianDistribution)
        @fact edge.distribution_type => GaussianDistribution
    end

    context("Edge constructor should throw an error if the FactorGraph is locked") do
        (node1, node2) = initializePairOfMockNodes()
        currentGraph().locked = true
        @fact_throws Edge(node1.i[:out], node2.i[:out])
    end

    context("Edge construction should couple interfaces to edge") do
        (node1, node2) = initializePairOfMockNodes()
        @fact node1.i[:out].edge => nothing
        @fact node2.i[:out].edge => nothing
        edge = Edge(node1.i[:out], node2.i[:out])
        @fact node1.i[:out].edge => edge
        @fact node2.i[:out].edge => edge
    end

    context("Edge should throw an error when the user attempts to reposition") do
        FactorGraph()
        node1 = TerminalNode(id=:node1)
        node2 = TerminalNode(id=:node2)
        node3 = TerminalNode(id=:node3)
        Edge(node1.i[:out], node2.i[:out])
        @fact_throws Edge(node1.i[:out], node3.i[:out])
    end

    context("Edges have a (pseudo) ordering") do
        FactorGraph()
        edge1 = Edge(TerminalNode().i[:out], TerminalNode().i[:out])
        edge2 = Edge(TerminalNode().i[:out], TerminalNode().i[:out])
        @fact (edge1 < edge2) => true
    end

    context("Edges have ids") do
        FactorGraph()
        TerminalNode(id=:a)
        TerminalNode(id=:b)
        my_edge = Edge(n(:a), n(:b), id="my_edge")
        @fact e(:my_edge) => my_edge

        TerminalNode(id=:c)
        TerminalNode(id=:d)
        my_edge2 = Edge(n(:c), n(:d))
        @fact my_edge2.id => :c_d
    end

    context("Edges can be sorted") do
        FactorGraph()
        node_a = TerminalNode(id=:a)
        node_b = FixedGainNode(id=:b)
        node_c = FixedGainNode(id=:c)
        node_d = TerminalNode(id=:d)
        edge_ab = Edge(node_a.i[:out], node_b.i[:in])
        edge_bc = Edge(node_b.i[:out], node_c.i[:in])
        edge_cd = Edge(node_c.i[:out], node_d.i[:out])
        sorted = sort!([edge_bc, edge_ab, edge_cd])
        @fact sorted => [edge_ab, edge_bc, edge_cd]
    end

    context("forwardMessage(), forwardMessage(), ensureMarginal!()") do
        FactorGraph()
        test_edge = Edge(MockNode(Message(GaussianDistribution())).i[:out], MockNode(Message(GaussianDistribution(m=3.0, V=2.0))).i[:out])
        @fact forwardMessage(test_edge) => Message(GaussianDistribution())
        @fact backwardMessage(test_edge) => Message(GaussianDistribution(m=3.0, V=2.0))
        test_edge.marginal = GaussianDistribution(m=3.0, V=2.0)
        @fact ForneyLab.ensureMarginal!(test_edge, GaussianDistribution) => GaussianDistribution(m=3.0, V=2.0)
        test_edge.marginal = nothing
        @fact ForneyLab.ensureMarginal!(test_edge, GaussianDistribution) => vague(GaussianDistribution)
    end
end

facts("Functions for collecting edges") do
    context("edges() should get all edges connected to the node set") do
        testnodes = initializeLoopyGraph()
        @fact edges(Set{Node}({testnodes[1], testnodes[2]})) => Set{Edge}({testnodes[1].i[:in].edge, testnodes[4].i[:in1].edge, testnodes[4].i[:out].edge})
    end
end
