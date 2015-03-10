#####################
# Unit tests
#####################

facts("Graph level unit tests") do
    context("FactorGraph() should initialize a factor graph") do
        fg = FactorGraph()
        @fact typeof(fg) => FactorGraph
        @fact length(fg.inference_schemes) => 1
        @fact typeof(fg.active_scheme) => InferenceScheme
        @fact fg.active_scheme => fg.inference_schemes[1]
        @fact fg.nodes => Set{Node}()
        @fact fg.edges => Set{Edge}()
        @fact current_graph => fg # Global should be set
    end

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

    context("InferenceScheme(::FactorGraph) should initialize an InferenceScheme and add it to the inference scheme array") do
        fg = FactorGraph()
        t1 = TerminalNode()
        t2 = TerminalNode()
        edge = Edge(t1.out, t2.out)
        scheme = InferenceScheme(fg) # Build new inference scheme and add to graph
        @fact typeof(scheme) => InferenceScheme
        @fact length(fg.inference_schemes) => 2 # Add to scheme list
        @fact fg.active_scheme => scheme # Set as active scheme
        @fact scheme.factorization[1].internal_edges => Set{Edge}({edge})
        @fact scheme.factorization[1].nodes => Set{Node}({t1, t2})
        @fact scheme.edge_to_subgraph => {edge => scheme.factorization[1]}
    end

    context("Subgraph() should initialize a subgraph") do
        scheme = InferenceScheme()
        sg = Subgraph()
        @fact typeof(sg) => Subgraph
        @fact typeof(sg.internal_schedule) => Schedule
        @fact typeof(sg.external_schedule) => ExternalSchedule
        @fact length(scheme.factorization) => 1 # Do not add to factorization
    end

    context("Subgraph(::InferenceScheme) should initialize a subgraph and add it to the inference scheme factorization") do
        scheme = InferenceScheme()
        sg = Subgraph(scheme)
        @fact typeof(sg) => Subgraph
        @fact length(scheme.factorization) => 2 # Add to factorization
    end

    context("currentGraph() should return a current graph object") do
        FactorGraph() # Reset
        @fact currentGraph() => current_graph
        my_graph = currentGraph()  # Make local pointer to global variable
        @fact typeof(my_graph) => FactorGraph
    end

    context("setCurrentGraph() should set a new current graph object") do
        my_first_graph = FactorGraph() # Reset
        my_second_graph = FactorGraph()
        @fact my_first_graph == current_graph => false
        setCurrentGraph(my_first_graph)
        @fact my_first_graph == current_graph => true
    end

    context("subgraph(edge) should return the subgraph where edge is internal") do
        (t1, a1, g1, t2, add1, g2) = initializeFactoringGraph()
        factorize!(Set{Edge}([t2.out.edge]))
        graph = currentGraph()
        @fact subgraph(t1.out.edge) => graph.active_scheme.factorization[1]
        @fact subgraph(t2.out.edge) => graph.active_scheme.factorization[2]
    end

    context("node(name::String) should return node with matching name") do
        graph_1 = FactorGraph()
        graph_2 = FactorGraph()
        n = MockNode(name="my_mocknode")
        Edge(n.out, MockNode().out)
        @fact_throws node("some_name")
        @fact_throws node("my_mocknode", graph_1)
        setCurrentGraph(graph_1)
        @fact_throws node("my_mocknode")
        @fact node("my_mocknode", graph_2) => n
        setCurrentGraph(graph_2)
        @fact node("my_mocknode") => n
    end
end

#####################
# Integration tests
#####################

facts("Graph level integration tests") do
    context("nodesConnectedToExternalEdges() should return all nodes (g) connected to external edges") do
        data = [1.0, 1.0, 1.0]

        # MF case
        (g_nodes, y_nodes, m_eq_nodes, gam_eq_nodes, q_m_edges, q_gam_edges, q_y_edges) = initializeGaussianNodeChain(data)
        n_sections = length(data)
        factorize!()
        graph = currentGraph()
        m_subgraph = subgraph(g_nodes[1].mean.edge)
        gam_subgraph = subgraph(g_nodes[1].precision.edge)
        y1_subgraph = subgraph(g_nodes[1].out.edge)
        y2_subgraph = subgraph(g_nodes[2].out.edge)
        y3_subgraph = subgraph(g_nodes[3].out.edge)
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
        m_gam_subgraph = subgraph(g_nodes[1].mean.edge)
        y1_subgraph = subgraph(g_nodes[1].out.edge)
        y2_subgraph = subgraph(g_nodes[2].out.edge)
        y3_subgraph = subgraph(g_nodes[3].out.edge)
        @fact Set(ForneyLab.nodesConnectedToExternalEdges(m_gam_subgraph)) => Set(g_nodes)
        @fact Set(ForneyLab.nodesConnectedToExternalEdges(y1_subgraph)) => Set([g_nodes[1]])
        @fact Set(ForneyLab.nodesConnectedToExternalEdges(y2_subgraph)) => Set([g_nodes[2]])
        @fact Set(ForneyLab.nodesConnectedToExternalEdges(y3_subgraph)) => Set([g_nodes[3]])
    end

    context("conformSubGraph!() should complete a subgraph with nodes and external edges based in its internal edges") do
        my_graph = FactorGraph()
        # On empty subgraph
        my_subgraph = my_graph.active_scheme.factorization[1]
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

    context("nodes() should return an array of all nodes in the graph") do
        # FactorGraph test
        testnodes = initializeLoopyGraph()
        found_nodes = nodes(currentGraph())
        @fact length(found_nodes) => length(testnodes)
        for node in testnodes
            @fact node in found_nodes => true
        end

        # Subgraph test
        found_nodes = nodes(currentGraph().active_scheme.factorization[1])
        @fact length(found_nodes) => length(testnodes)
        for node in testnodes
            @fact node in found_nodes => true
        end

        # Composite node test
        compnode = GainEqualityCompositeNode();
        found_nodes = nodes(compnode, depth=1)
        @fact length(found_nodes) => 2
        @fact compnode.equality_node in found_nodes => true
        @fact compnode.fixed_gain_node in found_nodes => true
    end

    context("edges() should get all edges internal (optionally external as well) to the argument node set") do
        testnodes = initializeLoopyGraph()
        @fact edges(Set{Node}({testnodes[1], testnodes[2]}), include_external=false) => Set{Edge}({testnodes[1].in1.edge})
        @fact edges(Set{Node}({testnodes[1], testnodes[2]})) => Set{Edge}({testnodes[1].in1.edge, testnodes[4].in1.edge, testnodes[4].out.edge})
    end

    context("edges() called on a subgraph should return all internal edges (optionally external as well) of the subgraph") do
        (t1, a1, g1, t2, t3) = initializeFactoringGraphWithoutLoop()
        graph = currentGraph()
        factorize!(graph)
        sg = graph.active_scheme.factorization[1]
        @fact edges(sg, include_external=false) => Set{Edge}({g1.variance.edge})
        @fact edges(sg) => Set{Edge}({g1.variance.edge, g1.out.edge, g1.mean.edge})
    end
end
