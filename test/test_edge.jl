#####################
# Integration tests
#####################

facts("Edge integration tests") do
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

    context("Edge constructor should write the distribution type") do
        (node1, node2) = initializePairOfMockNodes()
        edge = Edge(node1.out, node2.out, GaussianDistribution)
        @fact edge.distribution_type => GaussianDistribution
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
        FactorGraph()
        comp_node = GainEqualityCompositeNode()
        node = TerminalNode()
        edge = Edge(node.out, comp_node.in1)
        @fact comp_node.equality_node.interfaces[1].partner => node.out
        @fact comp_node.equality_node.interfaces[1].edge => edge
    end

    context("Edge should throw an error when the user attempts to reposition") do
        FactorGraph()
        node1 = TerminalNode(name="node1")
        node2 = TerminalNode(name="node2")
        node3 = TerminalNode(name="node3")
        Edge(node1.out, node2.out)
        @fact_throws Edge(node1.out, node3.out)
    end

    context("Edge should throw an error when the user attempts to extend upon a locked graph") do
        graph = FactorGraph()
        node1 = TerminalNode(name="node1")
        node2 = FixedGainNode(name="node2")
        node3 = TerminalNode(name="node3")
        Edge(node1.out, node2.in1)
        @fact graph.locked => false
        InferenceScheme() # Locks the graph structure (even though it is not done yet)
        @fact graph.locked => true
        @fact_throws Edge(node2.out, node3.out)
    end

    context("Edges have a (pseudo) ordering") do
        FactorGraph()
        edge1 = Edge(TerminalNode().out, TerminalNode().out)
        edge2 = Edge(TerminalNode().out, TerminalNode().out)
        @fact (edge1 < edge2) => true
    end

    context("Edges can be sorted") do
        FactorGraph()
        node_a = TerminalNode(name="a")
        node_b = FixedGainNode(name="b")
        node_c = FixedGainNode(name="c")
        node_d = TerminalNode(name="d")
        edge_ab = Edge(node_a.out, node_b.in1)
        edge_bc = Edge(node_b.out, node_c.in1)
        edge_cd = Edge(node_c.out, node_d.out)
        sorted = sort!([edge_bc, edge_ab, edge_cd])
        @fact sorted => [edge_ab, edge_bc, edge_cd]
    end

    context("forwardMessage(), forwardMessage(), ensureMarginal!()") do
        FactorGraph()
        test_edge = Edge(MockNode(Message(GaussianDistribution())).out, MockNode(Message(GaussianDistribution(m=3.0, V=2.0))).out)
        @fact forwardMessage(test_edge) => Message(GaussianDistribution())
        @fact backwardMessage(test_edge) => Message(GaussianDistribution(m=3.0, V=2.0))
        test_edge.marginal = GaussianDistribution(m=3.0, V=2.0)
        @fact ForneyLab.ensureMarginal!(test_edge, GaussianDistribution) => GaussianDistribution(m=3.0, V=2.0)
        test_edge.marginal = nothing
        @fact ForneyLab.ensureMarginal!(test_edge, GaussianDistribution) => vague(GaussianDistribution)
    end
end

facts("Functions for collecting edges") do
    context("edges() should get all edges internal (optionally external as well) to the argument node set") do
        testnodes = initializeLoopyGraph()
        @fact edges(Set{Node}({testnodes[1], testnodes[2]}), include_external=false) => Set{Edge}({testnodes[1].in1.edge})
        @fact edges(Set{Node}({testnodes[1], testnodes[2]})) => Set{Edge}({testnodes[1].in1.edge, testnodes[4].in1.edge, testnodes[4].out.edge})
    end

    context("edges() called on a subgraph should return all internal edges (optionally external as well) of the subgraph") do
        (t1, a1, g1, t2, t3) = initializeFactoringGraphWithoutLoop()
        scheme = InferenceScheme()
        factorize!(scheme)
        sg = scheme.factorization[1]
        @fact edges(sg, include_external=false) => Set{Edge}({g1.variance.edge})
        @fact edges(sg) => Set{Edge}({g1.variance.edge, g1.out.edge, g1.mean.edge})
    end
end
