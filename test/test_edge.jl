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

    context("Edge should throw an error when two interfaces of the same node are connected") do
        FactorGraph()
        node = GainNode()
        # Connect output directly to input
        @fact_throws Edge(node.interfaces[2], node.interfaces[1])
    end

    context("Edge should throw an error when a node is not in the current graph") do
        g1 = FactorGraph()
        node1 = MockNode()
        g2 = FactorGraph()
        node2 = MockNode()
        @fact_throws Edge(node1.i[:out], node2.i[:out])
    end

    context("Edge constructor should write the distribution type") do
        initializePairOfMockNodes()
        edge = Edge(n(:node1).i[:out], n(:node2).i[:out], GaussianDistribution)
        @fact edge.distribution_type --> GaussianDistribution
    end

    context("It is not possible to add an Edge to a locked FactorGraph or to deepcopy an Edge") do
        initializePairOfMockNodes()
        currentGraph().locked = true
        @fact_throws Edge(n(:node1).i[:out], n(:node2).i[:out])
        currentGraph().locked = false
        testedge = Edge(n(:node1).i[:out], n(:node2).i[:out])
        @fact_throws deepcopy(testedge)
    end

    context("Edge construction should couple partners and couple interfaces to edge") do
        initializePairOfMockNodes()
        @fact n(:node1).i[:out].edge --> nothing
        @fact n(:node2).i[:out].edge --> nothing
        @fact n(:node1).i[:out].partner --> nothing
        @fact n(:node2).i[:out].partner --> nothing
        edge = Edge(n(:node1).i[:out], n(:node2).i[:out])
        @fact n(:node1).i[:out].edge --> edge
        @fact n(:node2).i[:out].edge --> edge
        @fact n(:node1).i[:out].partner --> n(:node2).i[:out]
        @fact n(:node2).i[:out].partner --> n(:node1).i[:out]
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
        edge3 = Edge(TerminalNode().i[:out], TerminalNode().i[:out])
        edge4 = Edge(TerminalNode().i[:out], TerminalNode().i[:out])

        @fact (edge1 < edge2) --> true
        @fact ([edge1, edge3] < [edge2, edge4]) --> true
        @fact (edge1 < [edge2, edge4]) --> true
        @fact ([edge1, edge3] < edge2) --> true
    end

    context("Edges have ids") do
        FactorGraph()
        TerminalNode(id=:a)
        TerminalNode(id=:b)
        my_edge = Edge(n(:a), n(:b), id=:my_edge)
        @fact ForneyLab.e(:my_edge) --> my_edge

        TerminalNode(id=:c)
        TerminalNode(id=:d)
        my_edge2 = Edge(n(:c), n(:d))
        @fact my_edge2.id --> :c_d

        TerminalNode(id=:e)
        TerminalNode(id=:f)
        my_edge3 = Edge(n(:e), n(:f), id=:my_edge*3)
        @fact ForneyLab.e(:my_edge*3) --> my_edge3
    end

    context("Edges can be sorted") do
        FactorGraph()
        node_a = TerminalNode(id=:a)
        node_b = GainNode(id=:b)
        node_c = GainNode(id=:c)
        node_d = TerminalNode(id=:d)
        edge_ab = Edge(node_a.i[:out], node_b.i[:in])
        edge_bc = Edge(node_b.i[:out], node_c.i[:in])
        edge_cd = Edge(node_c.i[:out], node_d.i[:out])
        sorted = sort!([edge_bc, edge_ab, edge_cd])
        @fact sorted --> [edge_ab, edge_bc, edge_cd]
    end

    context("delete! should remove an edge and coupled write buffers") do
        g = initializePairOfMockNodes()
        edge = Edge(n(:node1).i[:out], n(:node2).i[:out])
        attachWriteBuffer(n(:node1).i[:out])
        attachWriteBuffer(ForneyLab.e(:node1_node2))

        delete!(g, edge)
        @fact length(g.edges) --> 0
        @fact n(:node1).i[:out].edge --> nothing
        @fact n(:node2).i[:out].edge --> nothing
        @fact n(:node1).i[:out].partner --> nothing
        @fact n(:node2).i[:out].partner --> nothing
        @fact length(g.write_buffers) --> 0
    end

    context("forwardMessage(), forwardMessage(), ensureMarginal!()") do
        FactorGraph()
        test_edge = Edge(MockNode().i[:out], MockNode().i[:out])
        test_edge.head.message = Message(GaussianDistribution(m=3.0, V=2.0))
        test_edge.tail.message = Message(GaussianDistribution())
        @fact forwardMessage(test_edge) --> Message(GaussianDistribution())
        @fact backwardMessage(test_edge) --> Message(GaussianDistribution(m=3.0, V=2.0))
        test_edge.marginal = GaussianDistribution(m=3.0, V=2.0)
        @fact ForneyLab.ensureMarginal!(test_edge, GaussianDistribution) --> GaussianDistribution(m=3.0, V=2.0)
        test_edge.marginal = nothing
        @fact ForneyLab.ensureMarginal!(test_edge, GaussianDistribution) --> vague(GaussianDistribution)
    end
end
