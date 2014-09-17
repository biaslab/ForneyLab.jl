#####################
# Unit tests
#####################

facts("Graph level unit tests") do
    context("FactorGraph() should initialize a factor graph") do
        fg = FactorGraph()
        @fact typeof(fg) => FactorGraph
        @fact length(fg.factorization) => 1
        @fact typeof(fg.factorization[1]) => Subgraph
        @fact typeof(fg.edge_to_subgraph) => Dict{Edge, Subgraph}
        @fact current_graph => fg # Global should be set
    end

    context("Subgraph() should initialize a subgraph and add it to the current graph") do
        sg = Subgraph()
        @fact typeof(sg) => Subgraph
        @fact typeof(sg.internal_schedule) => Schedule
        @fact typeof(sg.external_schedule) => ExternalSchedule
        graph = getCurrentGraph()
        @fact graph.factorization[2] => sg
    end

    context("getCurrentGraph() should return a current graph object") do
        FactorGraph() # Reset
        @fact getCurrentGraph() => current_graph
        my_graph = getCurrentGraph()  # Make local pointer to global variable
        @fact typeof(my_graph) => FactorGraph
    end

    context("setCurrentGraph() should set a new current graph object") do
        my_first_graph = FactorGraph() # Reset
        my_second_graph = FactorGraph()
        @fact my_first_graph == current_graph => false
        setCurrentGraph(my_first_graph)
        @fact my_first_graph == current_graph => true
    end

    context("getSubgraph(edge) should return the subgraph where edge is internal") do
        (driver, inhibitor, noise, add) = initializeLoopyGraph()
        factorize!(Set{Edge}({inhibitor.out.edge})) # Put this edge in a different subgraph
        graph = getCurrentGraph()
        @fact getSubgraph(driver.out.edge) => graph.factorization[1]
        @fact getSubgraph(inhibitor.out.edge) => graph.factorization[2]
    end
end

#####################
# Integration tests
#####################

facts("Graph level integration tests") do
    context("getNodesConnectedToExternalEdges() should return all nodes (g) connected to external edges") do
        data = [1.0, 1.0, 1.0]

        # MF case
        (g_nodes, y_nodes, m_eq_nodes, gam_eq_nodes, q_m_edges, q_gam_edges, q_y_edges) = initializeGaussianNodeChain(data)
        n_sections = length(data)
        factorizeMeanField!()
        graph = getCurrentGraph()
        m_subgraph = getSubgraph(g_nodes[1].mean.edge)
        gam_subgraph = getSubgraph(g_nodes[1].precision.edge)
        y1_subgraph = getSubgraph(g_nodes[1].out.edge)
        y2_subgraph = getSubgraph(g_nodes[2].out.edge)
        y3_subgraph = getSubgraph(g_nodes[3].out.edge)
        @fact Set(ForneyLab.getNodesConnectedToExternalEdges(m_subgraph)) => Set(g_nodes)
        @fact Set(ForneyLab.getNodesConnectedToExternalEdges(gam_subgraph)) => Set(g_nodes)
        @fact Set(ForneyLab.getNodesConnectedToExternalEdges(y1_subgraph)) => Set([g_nodes[1]])
        @fact Set(ForneyLab.getNodesConnectedToExternalEdges(y2_subgraph)) => Set([g_nodes[2]])
        @fact Set(ForneyLab.getNodesConnectedToExternalEdges(y3_subgraph)) => Set([g_nodes[3]])

        # Structured case
        (g_nodes, y_nodes, m_eq_nodes, gam_eq_nodes, q_m_edges, q_gam_edges, q_y_edges) = initializeGaussianNodeChain(data)
        n_sections = length(data)
        for edge in q_y_edges
            factorize!(Set{Edge}({edge}))
        end        
        graph = getCurrentGraph()
        m_gam_subgraph = getSubgraph(g_nodes[1].mean.edge)
        y1_subgraph = getSubgraph(g_nodes[1].out.edge)
        y2_subgraph = getSubgraph(g_nodes[2].out.edge)
        y3_subgraph = getSubgraph(g_nodes[3].out.edge)
        @fact Set(ForneyLab.getNodesConnectedToExternalEdges(m_gam_subgraph)) => Set(g_nodes)
        @fact Set(ForneyLab.getNodesConnectedToExternalEdges(y1_subgraph)) => Set([g_nodes[1]])
        @fact Set(ForneyLab.getNodesConnectedToExternalEdges(y2_subgraph)) => Set([g_nodes[2]])
        @fact Set(ForneyLab.getNodesConnectedToExternalEdges(y3_subgraph)) => Set([g_nodes[3]])
    end

    context("conformSubGraph!() should complete a subgraph with nodes and external edges based in its internal edges") do
        my_graph = FactorGraph()
        # On empty subgraph
        my_subgraph = my_graph.factorization[1]
        @fact length(my_subgraph.internal_edges) => 0
        ForneyLab.conformSubgraph!(my_subgraph)
        @fact length(my_subgraph.nodes) => 0
        @fact length(my_subgraph.external_edges) => 0
        # Initialize a subgraph
        node1 = MockNode()
        node2 = MockNode(2)
        node3 = MockNode() 
        edge1 = Edge(node1.out, node2.interfaces[1])
        edge2 = Edge(node2.interfaces[2], node3.out)
        @fact length(my_subgraph.internal_edges) => 2
        ForneyLab.conformSubgraph!(my_subgraph)
        @fact length(my_subgraph.nodes) => 3
        @fact length(my_subgraph.external_edges) => 0
        # Subgraph with external edges
        new_subgraph = Subgraph(Set{Node}(), Set{Edge}({edge2}), Set{Edge}(), Array(Interface, 0), Array(Node, 0))
        @fact length(new_subgraph.internal_edges) => 1
        ForneyLab.conformSubgraph!(new_subgraph)
        @fact length(new_subgraph.nodes) => 2
        @fact length(new_subgraph.external_edges) => 1
    end

    context("setUninformativeMarginals() should preset uninformative marginals at the appropriate places") do
        data = [1.0, 1.0, 1.0]

        # MF case
        (g_nodes, y_nodes, m_eq_nodes, gam_eq_nodes, q_m_edges, q_gam_edges, q_y_edges) = initializeGaussianNodeChain(data)
        n_sections = length(data)
        factorizeMeanField!()
        setUninformativeMarginals!()
        graph = getCurrentGraph()
        m_subgraph = getSubgraph(g_nodes[1].mean.edge)
        gam_subgraph = getSubgraph(g_nodes[1].precision.edge)
        y1_subgraph = getSubgraph(g_nodes[1].out.edge)
        y2_subgraph = getSubgraph(g_nodes[2].out.edge)
        y3_subgraph = getSubgraph(g_nodes[3].out.edge)

        @fact length(graph.approximate_marginals) => 9
        @fact graph.approximate_marginals[(g_nodes[1], m_subgraph)] => uninformative(GaussianDistribution)
        @fact graph.approximate_marginals[(g_nodes[2], m_subgraph)] => uninformative(GaussianDistribution)
        @fact graph.approximate_marginals[(g_nodes[3], m_subgraph)] => uninformative(GaussianDistribution)
        @fact graph.approximate_marginals[(g_nodes[1], gam_subgraph)] => uninformative(GammaDistribution)
        @fact graph.approximate_marginals[(g_nodes[2], gam_subgraph)] => uninformative(GammaDistribution)
        @fact graph.approximate_marginals[(g_nodes[3], gam_subgraph)] => uninformative(GammaDistribution)
        @fact graph.approximate_marginals[(g_nodes[1], y1_subgraph)] => 1.0
        @fact graph.approximate_marginals[(g_nodes[2], y2_subgraph)] => 1.0
        @fact graph.approximate_marginals[(g_nodes[3], y3_subgraph)] => 1.0

        # Structured case
        (g_nodes, y_nodes, m_eq_nodes, gam_eq_nodes, q_m_edges, q_gam_edges, q_y_edges) = initializeGaussianNodeChain(data)
        n_sections = length(data)
        for edge in q_y_edges
            factorize!(Set{Edge}({edge}))
        end        
        setUninformativeMarginals!()
        graph = getCurrentGraph()
        m_gam_subgraph = getSubgraph(g_nodes[1].mean.edge)
        y1_subgraph = getSubgraph(g_nodes[1].out.edge)
        y2_subgraph = getSubgraph(g_nodes[2].out.edge)
        y3_subgraph = getSubgraph(g_nodes[3].out.edge)

        @fact length(graph.approximate_marginals) => 6
        @fact graph.approximate_marginals[(g_nodes[1], m_gam_subgraph)] => uninformative(NormalGammaDistribution)
        @fact graph.approximate_marginals[(g_nodes[2], m_gam_subgraph)] => uninformative(NormalGammaDistribution)
        @fact graph.approximate_marginals[(g_nodes[3], m_gam_subgraph)] => uninformative(NormalGammaDistribution)
        @fact graph.approximate_marginals[(g_nodes[1], y1_subgraph)] => 1.0
        @fact graph.approximate_marginals[(g_nodes[2], y2_subgraph)] => 1.0
        @fact graph.approximate_marginals[(g_nodes[3], y3_subgraph)] => 1.0
    end

    context("getNodes() should return an array of all nodes in the graph") do
        nodes = initializeLoopyGraph()
        found_nodes = getNodes(getCurrentGraph())
        @fact length(found_nodes) => length(nodes) # FactorGraph test
        for node in nodes
            @fact node in found_nodes => true
        end

        found_nodes = getNodes(getCurrentGraph().factorization[1]) # Subgraph test
        @fact length(found_nodes) => length(nodes)
        for node in nodes
            @fact node in found_nodes => true
        end
    end

    context("getEdges() should get all edges internal (optionally external as well) to the argument node set") do
        nodes = initializeLoopyGraph()
        @fact getEdges(Set{Node}({nodes[1], nodes[2]}), include_external=false) => Set{Edge}({nodes[1].in1.edge})
        @fact getEdges(Set{Node}({nodes[1], nodes[2]})) => Set{Edge}({nodes[1].in1.edge, nodes[4].in1.edge, nodes[4].out.edge})
    end

    context("addChildNodes!() should add composite node's child nodes to the node array") do
        node = initializeGainEqualityCompositeNode(eye(1), false, [Message(GaussianDistribution()), Message(GaussianDistribution()), nothing])
        @fact ForneyLab.addChildNodes!(Set{Node}({node})) => Set{Node}({node, node.equality_node, node.fixed_gain_node})
    end
end