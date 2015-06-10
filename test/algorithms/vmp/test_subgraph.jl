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
        initializeFactoringGraphWithoutLoop()
        f = VMP.factorize()
        sg = f.factors[1]
        @fact nodes(sg) => Set{Node}({n(:g1), n(:t2)})
    end

    context("edges() called on a subgraph should return all internal edges (optionally external as well) of the subgraph") do
        initializeFactoringGraphWithoutLoop()
        f = VMP.factorize()
        sg = f.factors[1]
        @fact edges(sg, include_external=false) => Set{Edge}({n(:g1).i[:variance].edge})
        @fact edges(sg) => Set{Edge}({n(:g1).i[:variance].edge, n(:g1).i[:out].edge, n(:g1).i[:mean].edge})
    end
end

#####################
# Integration tests
#####################

facts("Subgraph integration tests") do
    context("externalEdges() should return all external edges") do
        data = [1.0, 1.0, 1.0]

        # MF case
        initializeGaussianNodeChain(data)
        n_sections = length(data)
        f = VMP.factorize()
        m_subgraph = f.edge_to_subgraph[n(:g1).i[:mean].edge]
        @fact VMP.externalEdges(m_subgraph) => Set({n(:g1).i[:out].edge, n(:g2).i[:out].edge, n(:g3).i[:out].edge, n(:g1).i[:precision].edge, n(:g2).i[:precision].edge, n(:g3).i[:precision].edge})

        # Structured case
        initializeGaussianNodeChain(data)
        n_sections = length(data)
        f = VMP.QFactorization()
        for edge in [e(:q_y1), e(:q_y2), e(:q_y3)]
            f = VMP.factorize!(Set{Edge}({edge}), f)
        end
        m_gam_subgraph = f.edge_to_subgraph[n(:g1).i[:mean].edge]
        @fact VMP.externalEdges(m_gam_subgraph) => Set({n(:g1).i[:out].edge, n(:g2).i[:out].edge, n(:g3).i[:out].edge})
    end

    context("nodesConnectedToExternalEdges() should return all nodes (g) connected to external edges") do
        data = [1.0, 1.0, 1.0]

        # MF case
        initializeGaussianNodeChain(data)
        n_sections = length(data)
        f = VMP.factorize()
        m_subgraph = f.edge_to_subgraph[n(:g1).i[:mean].edge]
        @fact VMP.nodesConnectedToExternalEdges(m_subgraph) => Set{Node}([n(:g1), n(:g2), n(:g3)])

        # Structured case
        initializeGaussianNodeChain(data)
        n_sections = length(data)
        f = VMP.QFactorization()
        for edge in [e(:q_y1), e(:q_y2), e(:q_y3)]
            f = VMP.factorize!(Set{Edge}({edge}), f)
        end
        m_gam_subgraph = f.edge_to_subgraph[n(:g1).i[:mean].edge]
        @fact VMP.nodesConnectedToExternalEdges(m_gam_subgraph) => Set{Node}([n(:g1), n(:g2), n(:g3)])
    end
end