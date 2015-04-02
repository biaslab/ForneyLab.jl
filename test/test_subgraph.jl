#####################
# Unit tests
#####################

facts("InferenceScheme unit tests") do
    context("Subgraph() should initialize a subgraph") do
        scheme = InferenceScheme()
        sg = Subgraph()
        @fact typeof(sg) => Subgraph
        @fact typeof(sg.internal_schedule) => Schedule
        @fact typeof(sg.external_schedule) => ExternalSchedule
        @fact length(scheme.factorization) => 1 # Do not add to factorization
    end

    context("subgraph(edge) should return the subgraph where edge is internal") do
        (t1, a1, g1, t2, add1, g2) = initializeFactoringGraph()
        scheme = InferenceScheme()
        factorize!(Set{Edge}([t2.out.edge]))
        graph = currentGraph()
        @fact subgraph(scheme, t1.out.edge) => scheme.factorization[1]
        @fact subgraph(scheme, t2.out.edge) => scheme.factorization[2]
    end
end


#####################
# Integration tests
#####################

facts("Subgraph integration tests") do
    context("nodesConnectedToExternalEdges() should return all nodes (g) connected to external edges") do
        data = [1.0, 1.0, 1.0]

        # MF case
        (g_nodes, y_nodes, m_eq_nodes, gam_eq_nodes, q_m_edges, q_gam_edges, q_y_edges) = initializeGaussianNodeChain(data)
        n_sections = length(data)
        factorize!()
        graph = currentGraph()
        scheme = currentScheme()
        m_subgraph = subgraph(scheme, g_nodes[1].mean.edge)
        gam_subgraph = subgraph(scheme, g_nodes[1].precision.edge)
        y1_subgraph = subgraph(scheme, g_nodes[1].out.edge)
        y2_subgraph = subgraph(scheme, g_nodes[2].out.edge)
        y3_subgraph = subgraph(scheme, g_nodes[3].out.edge)
        @fact Set(ForneyLab.nodesConnectedToExternalEdges(m_subgraph)) => Set(g_nodes)
        @fact Set(ForneyLab.nodesConnectedToExternalEdges(gam_subgraph)) => Set(g_nodes)
        @fact Set(ForneyLab.nodesConnectedToExternalEdges(y1_subgraph)) => Set([g_nodes[1]])
        @fact Set(ForneyLab.nodesConnectedToExternalEdges(y2_subgraph)) => Set([g_nodes[2]])
        @fact Set(ForneyLab.nodesConnectedToExternalEdges(y3_subgraph)) => Set([g_nodes[3]])

        # Structured case
        (g_nodes, y_nodes, m_eq_nodes, gam_eq_nodes, q_m_edges, q_gam_edges, q_y_edges) = initializeGaussianNodeChain(data)
        n_sections = length(data)
        for edge in q_y_edges
            factorize!(Set{Edge}({edge}))
        end
        graph = currentGraph()
        scheme = currentScheme()
        m_gam_subgraph = subgraph(scheme, g_nodes[1].mean.edge)
        y1_subgraph = subgraph(scheme, g_nodes[1].out.edge)
        y2_subgraph = subgraph(scheme, g_nodes[2].out.edge)
        y3_subgraph = subgraph(scheme, g_nodes[3].out.edge)
        @fact Set(ForneyLab.nodesConnectedToExternalEdges(m_gam_subgraph)) => Set(g_nodes)
        @fact Set(ForneyLab.nodesConnectedToExternalEdges(y1_subgraph)) => Set([g_nodes[1]])
        @fact Set(ForneyLab.nodesConnectedToExternalEdges(y2_subgraph)) => Set([g_nodes[2]])
        @fact Set(ForneyLab.nodesConnectedToExternalEdges(y3_subgraph)) => Set([g_nodes[3]])
    end

    context("conformSubGraph!() should complete a subgraph with nodes and external edges based in its internal edges") do
        my_graph = FactorGraph()
        node1 = MockNode()
        node2 = MockNode(2)
        node3 = MockNode()
        edge1 = Edge(node1.out, node2.interfaces[1])
        edge2 = Edge(node2.interfaces[2], node3.out)

        scheme = InferenceScheme()
        my_subgraph = scheme.factorization[1]
        ForneyLab.conformSubgraph!(my_subgraph)
        @fact length(my_subgraph.internal_edges) => 2
        @fact length(my_subgraph.nodes) => 3
        @fact length(my_subgraph.external_edges) => 0

        # Subgraph with external edges
        new_subgraph = Subgraph(Set{Node}(), Set{Edge}({edge2}), Set{Edge}(), Array(Interface, 0), Array(Node, 0))
        ForneyLab.conformSubgraph!(new_subgraph)
        @fact length(new_subgraph.internal_edges) => 1
        @fact length(new_subgraph.nodes) => 2
        @fact length(new_subgraph.external_edges) => 1
    end
end