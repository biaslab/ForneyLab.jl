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

    context("It is not possible to deepcopy an Edge") do
        initializePairOfMockNodes()
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
        @fact eg(:my_edge) --> my_edge

        TerminalNode(id=:c)
        TerminalNode(id=:d)
        my_edge2 = Edge(n(:c), n(:d))
        @fact my_edge2.id --> :c_d

        TerminalNode(id=:e)
        TerminalNode(id=:f)
        my_edge3 = Edge(n(:e), n(:f), id=:my_edge*3)
        @fact eg(:my_edge*3) --> my_edge3
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

    context("forwardMessage(), forwardMessage(), ensureMarginal!()") do
        FactorGraph()
        test_edge = Edge(MockNode().i[:out], MockNode().i[:out])
        test_edge.head.message = Message(Gaussian(m=3.0, V=2.0))
        test_edge.tail.message = Message(Gaussian())
        @fact forwardMessage(test_edge) --> Message(Gaussian())
        @fact backwardMessage(test_edge) --> Message(Gaussian(m=3.0, V=2.0))
        test_edge.marginal = Gaussian(m=3.0, V=2.0)
        @fact ForneyLab.ensureMarginal!(test_edge, Gaussian) --> Gaussian(m=3.0, V=2.0)
        test_edge.marginal = nothing
        @fact ForneyLab.ensureMarginal!(test_edge, Gaussian) --> vague(Gaussian)
    end

    context("calculateMarginal, calculateMarginal!") do
        FactorGraph()
        test_edge = Edge(MockNode().i[:out], MockNode().i[:out])
        @fact_throws calculateMarginal(test_edge)
        @fact_throws calculateMarginal!(test_edge)
        d1 = Gaussian(m=3.0, V=1.0)
        d2 = Gaussian()
        marg = d1 * d2
        test_edge.head.message = Message(d1)
        test_edge.tail.message = Message(d2)
        @fact calculateMarginal(test_edge) --> marg
        @fact test_edge.marginal --> nothing
        @fact calculateMarginal!(test_edge) --> marg
        @fact test_edge.marginal --> marg
    end
end
