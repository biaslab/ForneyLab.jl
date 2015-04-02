#####################
# Unit tests
#####################

facts("InferenceScheme unit tests") do
    context("Subgraph() should initialize a subgraph") do
        sg = Subgraph()
        @fact typeof(sg) => Subgraph
        @fact typeof(sg.internal_edges) => Set{Edge}
        @fact typeof(sg.internal_schedule) => Schedule
    end

    context("subgraph(edge) should return the subgraph where edge is internal") do
        (t1, a1, g1, t2, add1, g2) = initializeFactoringGraph()
        graph = currentGraph()
        scheme = InferenceScheme()
        factorize!(Set{Edge}([t2.out.edge]))
        @fact subgraph(scheme, t1.out.edge) => scheme.factorization[1]
        @fact subgraph(scheme, t2.out.edge) => scheme.factorization[2]
    end
end


#####################
# Integration tests
#####################

facts("Subgraph integration tests") do
    # context("externalEdges() should return all external edges") do
    #     data = [1.0, 1.0, 1.0]

    #     # MF case
    #     (g_nodes, y_nodes, m_eq_nodes, gam_eq_nodes, q_m_edges, q_gam_edges, q_y_edges) = initializeGaussianNodeChain(data)
    #     n_sections = length(data)
    #     factorize!()
    #     graph = currentGraph()
    #     scheme = currentScheme()
    #     m_subgraph = subgraph(scheme, g_nodes[1].mean.edge)
    #     gam_subgraph = subgraph(scheme, g_nodes[1].precision.edge)
    #     y1_subgraph = subgraph(scheme, g_nodes[1].out.edge)
    #     y2_subgraph = subgraph(scheme, g_nodes[2].out.edge)
    #     y3_subgraph = subgraph(scheme, g_nodes[3].out.edge)
    #     @fact ForneyLab.externalEdges(m_subgraph) => Set(g_nodes)
    #     @fact ForneyLab.externalEdges(gam_subgraph) => Set(g_nodes)
    #     @fact ForneyLab.externalEdges(y1_subgraph) => Set([g_nodes[1]])
    #     @fact ForneyLab.externalEdges(y2_subgraph) => Set([g_nodes[2]])
    #     @fact ForneyLab.externalEdges(y3_subgraph) => Set([g_nodes[3]])

    #     # Structured case
    #     (g_nodes, y_nodes, m_eq_nodes, gam_eq_nodes, q_m_edges, q_gam_edges, q_y_edges) = initializeGaussianNodeChain(data)
    #     n_sections = length(data)
    #     for edge in q_y_edges
    #         factorize!(Set{Edge}({edge}))
    #     end
    #     graph = currentGraph()
    #     scheme = currentScheme()
    #     m_gam_subgraph = subgraph(scheme, g_nodes[1].mean.edge)
    #     y1_subgraph = subgraph(scheme, g_nodes[1].out.edge)
    #     y2_subgraph = subgraph(scheme, g_nodes[2].out.edge)
    #     y3_subgraph = subgraph(scheme, g_nodes[3].out.edge)
    #     @fact ForneyLab.externalEdges(m_gam_subgraph) => Set(g_nodes)
    #     @fact ForneyLab.externalEdges(y1_subgraph) => Set([g_nodes[1]])
    #     @fact ForneyLab.externalEdges(y2_subgraph) => Set([g_nodes[2]])
    #     @fact ForneyLab.externalEdges(y3_subgraph) => Set([g_nodes[3]])
    # end

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
        @fact ForneyLab.nodesConnectedToExternalEdges(m_subgraph) => Set(g_nodes)
        @fact ForneyLab.nodesConnectedToExternalEdges(gam_subgraph) => Set(g_nodes)
        @fact ForneyLab.nodesConnectedToExternalEdges(y1_subgraph) => Set([g_nodes[1]])
        @fact ForneyLab.nodesConnectedToExternalEdges(y2_subgraph) => Set([g_nodes[2]])
        @fact ForneyLab.nodesConnectedToExternalEdges(y3_subgraph) => Set([g_nodes[3]])

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
        @fact ForneyLab.nodesConnectedToExternalEdges(m_gam_subgraph) => Set(g_nodes)
        @fact ForneyLab.nodesConnectedToExternalEdges(y1_subgraph) => Set([g_nodes[1]])
        @fact ForneyLab.nodesConnectedToExternalEdges(y2_subgraph) => Set([g_nodes[2]])
        @fact ForneyLab.nodesConnectedToExternalEdges(y3_subgraph) => Set([g_nodes[3]])
    end
end