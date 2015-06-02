#####################
# Integration tests
#####################

facts("Edge integration tests") do
    context("Nodes can be coupled by edges using the explicit interface names") do
        initializePairOfNodes()
        # Couple the interfaces that carry GeneralMessage
        edge = Edge(n(:node2).i[:out], n(:node1).i[:in]) # Edge from node 2 to node 1
        testInterfaceConnections(n(:node1), n(:node2))
    end

    context("Edge should throw an error when two interfaces on the same node are connected") do
        node = FixedGainNode()
        # Connect output directly to input
        @fact_throws Edge(node.interfaces[2], node.interfaces[1])
    end

    context("Edge constructor should write the distribution type") do
        initializePairOfMockNodes()
        edge = Edge(n(:node1).i[:out], n(:node2).i[:out], GaussianDistribution)
        @fact edge.distribution_type => GaussianDistribution
    end

    context("Edge constructor should throw an error if the FactorGraph is locked") do
        initializePairOfMockNodes()
        currentGraph().locked = true
        @fact_throws Edge(n(:node1).i[:out], n(:node2).i[:out])
    end

    context("Edge construction should couple interfaces to edge") do
        initializePairOfMockNodes()
        @fact n(:node1).i[:out].edge => nothing
        @fact n(:node2).i[:out].edge => nothing
        edge = Edge(n(:node1).i[:out], n(:node2).i[:out])
        @fact n(:node1).i[:out].edge => edge
        @fact n(:node2).i[:out].edge => edge
    end

    context("Edge should throw an error when the user attempts to reposition") do
        FactorGraph()
        TerminalNode(id=:node1)
        TerminalNode(id=:node2)
        TerminalNode(id=:node3)
        Edge(n(:node1).i[:out], n(:node2).i[:out])
        @fact_throws Edge(n(:node1).i[:out], n(:node3).i[:out])
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
        my_edge = Edge(n(:a), n(:b), id=:my_edge)
        @fact e(:my_edge) => my_edge

        TerminalNode(id=:c)
        TerminalNode(id=:d)
        my_edge2 = Edge(n(:c), n(:d))
        @fact my_edge2.id => :c_d

        TerminalNode(id=:e)
        TerminalNode(id=:f)
        my_edge3 = Edge(n(:e), n(:f), id=s(:my_edge,3))
        @fact e(:my_edge,3) => my_edge3
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
        initializeLoopyGraph()
        @fact edges(Set{Node}({n(:driver), n(:inhibitor)})) => Set{Edge}({n(:driver).i[:in].edge, n(:add).i[:in1].edge, n(:add).i[:out].edge})
    end
end
