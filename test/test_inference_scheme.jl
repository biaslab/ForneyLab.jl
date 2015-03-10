#####################
# Unit tests
#####################

facts("InferenceScheme unit tests") do
    context("InferenceScheme() should initialize an InferenceScheme") do
        fg = FactorGraph()
        scheme = InferenceScheme()
        @fact typeof(scheme) => InferenceScheme
        @fact length(fg.inference_schemes) => 1 # Do not add to scheme list
        @fact fg.active_scheme => fg.inference_schemes[1] # Do not set as active scheme
        @fact typeof(scheme.factorization[1]) => Subgraph
        @fact typeof(scheme.edge_to_subgraph) => Dict{Edge, Subgraph}
        @fact typeof(scheme.approximate_marginals) => Dict{(Node, Subgraph), ProbabilityDistribution}
        @fact typeof(scheme.read_buffers) => Dict{TerminalNode, Vector}
        @fact typeof(scheme.write_buffers) => Dict{Union(Edge,Interface), Vector}
        @fact typeof(scheme.time_wraps) => Vector{(TerminalNode, TerminalNode)}
    end

    context("InferenceScheme!(::FactorGraph) should initialize an InferenceScheme and add it to the inference scheme array") do
        fg = FactorGraph()
        t1 = TerminalNode()
        t2 = TerminalNode()
        edge = Edge(t1.out, t2.out)
        scheme = InferenceScheme!(fg) # Build new inference scheme and add to graph
        @fact typeof(scheme) => InferenceScheme
        @fact length(fg.inference_schemes) => 2 # Add to scheme list
        @fact fg.active_scheme => scheme # Set as active scheme
        @fact scheme.factorization[1].internal_edges => Set{Edge}({edge})
        @fact scheme.factorization[1].nodes => Set{Node}({t1, t2})
        @fact scheme.edge_to_subgraph => {edge => scheme.factorization[1]}
    end
end


#####################
# Integration tests
#####################

facts("InferenceScheme integration tests") do
    context("setVagueMarginals() should set vague marginals at the appropriate places") do
        data = [1.0, 1.0, 1.0]

        # MF case
        (g_nodes, y_nodes, m_eq_nodes, gam_eq_nodes, q_m_edges, q_gam_edges, q_y_edges) = initializeGaussianNodeChain(data)
        n_sections = length(data)
        factorize!()
        setVagueMarginals!()
        graph = currentGraph()
        m_subgraph = subgraph(g_nodes[1].mean.edge)
        gam_subgraph = subgraph(g_nodes[1].precision.edge)
        y1_subgraph = subgraph(g_nodes[1].out.edge)
        y2_subgraph = subgraph(g_nodes[2].out.edge)
        y3_subgraph = subgraph(g_nodes[3].out.edge)

        @fact length(graph.active_scheme.approximate_marginals) => 9
        @fact graph.active_scheme.approximate_marginals[(g_nodes[1], m_subgraph)] => vague(GaussianDistribution)
        @fact graph.active_scheme.approximate_marginals[(g_nodes[2], m_subgraph)] => vague(GaussianDistribution)
        @fact graph.active_scheme.approximate_marginals[(g_nodes[3], m_subgraph)] => vague(GaussianDistribution)
        @fact graph.active_scheme.approximate_marginals[(g_nodes[1], gam_subgraph)] => vague(GammaDistribution)
        @fact graph.active_scheme.approximate_marginals[(g_nodes[2], gam_subgraph)] => vague(GammaDistribution)
        @fact graph.active_scheme.approximate_marginals[(g_nodes[3], gam_subgraph)] => vague(GammaDistribution)
        @fact graph.active_scheme.approximate_marginals[(g_nodes[1], y1_subgraph)] => vague(GaussianDistribution)
        @fact graph.active_scheme.approximate_marginals[(g_nodes[2], y2_subgraph)] => vague(GaussianDistribution)
        @fact graph.active_scheme.approximate_marginals[(g_nodes[3], y3_subgraph)] => vague(GaussianDistribution)

        # Structured case
        (g_nodes, y_nodes, m_eq_nodes, gam_eq_nodes, q_m_edges, q_gam_edges, q_y_edges) = initializeGaussianNodeChain(data)
        n_sections = length(data)
        for edge in q_y_edges
            factorize!(Set{Edge}({edge}))
        end
        setVagueMarginals!()
        graph = currentGraph()
        m_gam_subgraph = subgraph(g_nodes[1].mean.edge)
        y1_subgraph = subgraph(g_nodes[1].out.edge)
        y2_subgraph = subgraph(g_nodes[2].out.edge)
        y3_subgraph = subgraph(g_nodes[3].out.edge)

        @fact length(graph.active_scheme.approximate_marginals) => 6
        @fact graph.active_scheme.approximate_marginals[(g_nodes[1], m_gam_subgraph)] => vague(NormalGammaDistribution)
        @fact graph.active_scheme.approximate_marginals[(g_nodes[2], m_gam_subgraph)] => vague(NormalGammaDistribution)
        @fact graph.active_scheme.approximate_marginals[(g_nodes[3], m_gam_subgraph)] => vague(NormalGammaDistribution)
        @fact graph.active_scheme.approximate_marginals[(g_nodes[1], y1_subgraph)] => vague(GaussianDistribution)
        @fact graph.active_scheme.approximate_marginals[(g_nodes[2], y2_subgraph)] => vague(GaussianDistribution)
        @fact graph.active_scheme.approximate_marginals[(g_nodes[3], y3_subgraph)] => vague(GaussianDistribution)
    end
end