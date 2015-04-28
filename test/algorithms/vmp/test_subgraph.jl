#####################
# Unit tests
#####################

facts("Subgraph unit tests") do
    context("Subgraph() should initialize a subgraph") do
        sg = VMP.Subgraph()
        @fact typeof(sg) => VMP.Subgraph
        @fact typeof(sg.internal_edges) => Set{Edge}
        @fact typeof(sg.internal_schedule) => Schedule
        @fact typeof(sg.external_schedule) => Array{Node, 1}
    end
end

facts("Nodes and edges overloadings for Subgraph") do
    context("nodes() called on a subgraph should return all internal nodes of the subgraph") do
        (t1, a1, g1, t2, t3) = initializeFactoringGraphWithoutLoop()
        f = VMP.factorize()
        sg = f.factors[1]
        @fact nodes(sg) => Set{Node}({g1, t2})
    end

    context("edges() called on a subgraph should return all internal edges (optionally external as well) of the subgraph") do
        (t1, a1, g1, t2, t3) = initializeFactoringGraphWithoutLoop()
        f = VMP.factorize()
        sg = f.factors[1]
        @fact edges(sg, include_external=false) => Set{Edge}({g1.variance.edge})
        @fact edges(sg) => Set{Edge}({g1.variance.edge, g1.out.edge, g1.mean.edge})
    end
end

#####################
# Integration tests
#####################

facts("Subgraph integration tests") do
    context("externalEdges() should return all external edges") do
        data = [1.0, 1.0, 1.0]

        # MF case
        (g_nodes, y_nodes, m_eq_nodes, gam_eq_nodes, q_m_edges, q_gam_edges, q_y_edges) = initializeGaussianNodeChain(data)
        n_sections = length(data)
        f = VMP.factorize()
        m_subgraph = f.edge_to_subgraph[g_nodes[1].mean.edge]
        @fact VMP.externalEdges(m_subgraph) => Set({g_nodes[1].out.edge, g_nodes[2].out.edge, g_nodes[3].out.edge, g_nodes[1].precision.edge, g_nodes[2].precision.edge, g_nodes[3].precision.edge})

        # Structured case
        (g_nodes, y_nodes, m_eq_nodes, gam_eq_nodes, q_m_edges, q_gam_edges, q_y_edges) = initializeGaussianNodeChain(data)
        n_sections = length(data)
        f = VMP.QFactorization()
        for edge in q_y_edges
            f = VMP.factorize!(Set{Edge}({edge}), f)
        end
        m_gam_subgraph = f.edge_to_subgraph[g_nodes[1].mean.edge]
        @fact VMP.externalEdges(m_gam_subgraph) => Set({g_nodes[1].out.edge, g_nodes[2].out.edge, g_nodes[3].out.edge})
    end

    context("nodesConnectedToExternalEdges() should return all nodes (g) connected to external edges") do
        data = [1.0, 1.0, 1.0]

        # MF case
        (g_nodes, y_nodes, m_eq_nodes, gam_eq_nodes, q_m_edges, q_gam_edges, q_y_edges) = initializeGaussianNodeChain(data)
        n_sections = length(data)
        f = VMP.factorize()
        m_subgraph = f.edge_to_subgraph[g_nodes[1].mean.edge]
        @fact VMP.nodesConnectedToExternalEdges(m_subgraph) => Set(g_nodes)

        # Structured case
        (g_nodes, y_nodes, m_eq_nodes, gam_eq_nodes, q_m_edges, q_gam_edges, q_y_edges) = initializeGaussianNodeChain(data)
        n_sections = length(data)
        f = VMP.QFactorization()
        for edge in q_y_edges
            f = VMP.factorize!(Set{Edge}({edge}), f)
        end
        m_gam_subgraph = f.edge_to_subgraph[g_nodes[1].mean.edge]
        @fact VMP.nodesConnectedToExternalEdges(m_gam_subgraph) => Set(g_nodes)
    end
end