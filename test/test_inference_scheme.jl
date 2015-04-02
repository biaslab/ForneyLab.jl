#####################
# Unit tests
#####################

facts("InferenceScheme unit tests") do
    context("InferenceScheme() should initialize an InferenceScheme") do
        fg = FactorGraph()
        scheme = InferenceScheme()
        @fact typeof(scheme) => InferenceScheme
        @fact currentScheme() => scheme
        @fact typeof(scheme.factorization[1]) => Subgraph
        @fact typeof(scheme.edge_to_subgraph) => Dict{Edge, Subgraph}
        @fact typeof(scheme.q_distributions) => Dict{Set{Edge}, ProbabilityDistribution}
        @fact typeof(scheme.read_buffers) => Dict{TerminalNode, Vector}
        @fact typeof(scheme.write_buffers) => Dict{Union(Edge,Interface), Vector}
        @fact typeof(scheme.time_wraps) => Vector{(TerminalNode, TerminalNode)}
    end
end


#####################
# Integration tests
#####################

facts("InferenceScheme integration tests") do
    context("setVagueQDistributions() should set vague marginals at the appropriate places") do
        data = [1.0, 1.0, 1.0]

        # MF case
        (g_nodes, y_nodes, m_eq_nodes, gam_eq_nodes, q_m_edges, q_gam_edges, q_y_edges) = initializeGaussianNodeChain(data)
        n_sections = length(data)
        scheme = InferenceScheme()
        factorize!()
        setVagueQDistributions!()
        graph = currentGraph()
        m_subgraph = subgraph(scheme, g_nodes[1].mean.edge)
        gam_subgraph = subgraph(scheme, g_nodes[1].precision.edge)
        y1_subgraph = subgraph(scheme, g_nodes[1].out.edge)
        y2_subgraph = subgraph(scheme, g_nodes[2].out.edge)
        y3_subgraph = subgraph(scheme, g_nodes[3].out.edge)

        @fact length(scheme.q_distributions) => 9
        @fact qDistribution(scheme, g_nodes[1], m_subgraph) => vague(GaussianDistribution)
        @fact qDistribution(scheme, g_nodes[2], m_subgraph) => vague(GaussianDistribution)
        @fact qDistribution(scheme, g_nodes[3], m_subgraph) => vague(GaussianDistribution)
        @fact qDistribution(scheme, g_nodes[1], gam_subgraph) => vague(GammaDistribution)
        @fact qDistribution(scheme, g_nodes[2], gam_subgraph) => vague(GammaDistribution)
        @fact qDistribution(scheme, g_nodes[3], gam_subgraph) => vague(GammaDistribution)
        @fact qDistribution(scheme, g_nodes[1], y1_subgraph) => vague(GaussianDistribution)
        @fact qDistribution(scheme, g_nodes[2], y2_subgraph) => vague(GaussianDistribution)
        @fact qDistribution(scheme, g_nodes[3], y3_subgraph) => vague(GaussianDistribution)

        # Structured case
        (g_nodes, y_nodes, m_eq_nodes, gam_eq_nodes, q_m_edges, q_gam_edges, q_y_edges) = initializeGaussianNodeChain(data)
        n_sections = length(data)
        scheme = InferenceScheme()
        for edge in q_y_edges
            factorize!(Set{Edge}({edge}))
        end
        setVagueQDistributions!()
        graph = currentGraph()
        m_gam_subgraph = subgraph(scheme, g_nodes[1].mean.edge)
        y1_subgraph = subgraph(scheme, g_nodes[1].out.edge)
        y2_subgraph = subgraph(scheme, g_nodes[2].out.edge)
        y3_subgraph = subgraph(scheme, g_nodes[3].out.edge)

        @fact length(scheme.q_distributions) => 6
        @fact qDistribution(scheme, g_nodes[1], m_gam_subgraph) => vague(NormalGammaDistribution)
        @fact qDistribution(scheme, g_nodes[2], m_gam_subgraph) => vague(NormalGammaDistribution)
        @fact qDistribution(scheme, g_nodes[3], m_gam_subgraph) => vague(NormalGammaDistribution)
        @fact qDistribution(scheme, g_nodes[1], y1_subgraph) => vague(GaussianDistribution)
        @fact qDistribution(scheme, g_nodes[2], y2_subgraph) => vague(GaussianDistribution)
        @fact qDistribution(scheme, g_nodes[3], y3_subgraph) => vague(GaussianDistribution)
    end
end