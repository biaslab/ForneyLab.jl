#####################
# Integration tests
#####################

facts("QFactorization integration tests") do
    context("extend() should extend a set of edges to envelope deterministic nodes") do
        (t1, a1, g1, t2, add1, g2) = initializeFactoringGraph()
        cluster = VMP.extend(Set{Edge}([g1.i[:out].edge, g1.i[:mean].edge]))
        @fact cluster => Set{Edge}({t1.i[:out].edge, a1.i[:out].edge, g1.i[:out].edge, add1.i[:out].edge, g2.i[:out].edge})
    end

    context("factorize!()") do
        context("Should include argument edges in a new subgraph") do
            (t1, a1, g1, t2, add1, g2) = initializeFactoringGraph()
            f = VMP.factorize!(Set{Edge}([t2.i[:out].edge]))
            @fact f.factors[2].internal_edges => Set{Edge}({t2.i[:out].edge})
            @fact f.factors[1].internal_edges => Set{Edge}({t1.i[:out].edge, a1.i[:out].edge, g1.i[:out].edge, add1.i[:out].edge, g2.i[:out].edge})
        end

        context("Should update the edge_to_subgraph mapping for the graph") do
            (t1, a1, g1, t2, add1, g2) = initializeFactoringGraph()
            f = VMP.factorize!(Set{Edge}([t2.i[:out].edge]))
            @fact f.edge_to_subgraph[t1.i[:out].edge] => f.factors[1]
            @fact f.edge_to_subgraph[a1.i[:out].edge] => f.factors[1]
            @fact f.edge_to_subgraph[g1.i[:out].edge] => f.factors[1]
            @fact f.edge_to_subgraph[t2.i[:out].edge] => f.factors[2]
            @fact f.edge_to_subgraph[add1.i[:out].edge] => f.factors[1]
            @fact f.edge_to_subgraph[g2.i[:out].edge] => f.factors[1]
        end

        context("Should not factorize a GaussianNode with fixed mean and variance") do
            (t, gain, gauss) = initializeGaussianFactoringGraph()
            f = VMP.factorize!(Set{Edge}([t.i[:out].edge]))
            @fact length(f.factors) => 1
            @fact f.edge_to_subgraph[t.i[:out].edge] => f.factors[1]
            @fact f.edge_to_subgraph[gauss.i[:out].edge] => f.factors[1]
        end

        context("Should retain node names") do
            (t1, gain, t2) = initializeSimpleFactoringGraph()
            f = VMP.factorize!(Set{Edge}([t1.i[:out].edge]))
            @fact node("t1") => t1
        end
    end

    context("factorize() should output a mean field factorized graph") do
        data = [1.0, 1.0, 1.0]
        (g_nodes, y_nodes, m_eq_nodes, gam_eq_nodes, q_m_edges, q_gam_edges, q_y_edges) = initializeGaussianNodeChain(data)
        f = VMP.factorize()
        gam_set = Set{Edge}()
        for gam_eq_node in gam_eq_nodes
            for interface in gam_eq_node.interfaces
                push!(gam_set, interface.edge)
            end
        end
        m_set = Set{Edge}()
        for m_eq_node in m_eq_nodes
            for interface in m_eq_node.interfaces
                push!(m_set, interface.edge)
            end
        end
        @fact length(f.factors) => 5

        @fact f.edge_to_subgraph[g_nodes[1].i[:mean].edge].internal_edges => m_set 
        @fact f.edge_to_subgraph[g_nodes[1].i[:precision].edge].internal_edges => gam_set 
        @fact f.edge_to_subgraph[g_nodes[1].i[:out].edge].internal_edges => Set{Edge}([q_y_edges[1]]) 
        @fact f.edge_to_subgraph[g_nodes[2].i[:out].edge].internal_edges => Set{Edge}([q_y_edges[2]])
        @fact f.edge_to_subgraph[g_nodes[3].i[:out].edge].internal_edges => Set{Edge}([q_y_edges[3]])
    end
end

facts("vagueQDistributions() should set vague marginals at the appropriate places") do
    data = [1.0, 1.0, 1.0]

    # MF case
    (g_nodes, y_nodes, m_eq_nodes, gam_eq_nodes, q_m_edges, q_gam_edges, q_y_edges) = initializeGaussianNodeChain(data)
    n_sections = length(data)

    f = VMP.factorize()
    VMP.generateSchedule!(f) # Generate and store internal and external schedules on factorization subgraphs
    qs = VMP.vagueQDistributions(f) # Initialize vague q distributions

    m_subgraph = f.edge_to_subgraph[g_nodes[1].i[:mean].edge]
    gam_subgraph = f.edge_to_subgraph[g_nodes[1].i[:precision].edge]
    y1_subgraph = f.edge_to_subgraph[g_nodes[1].i[:out].edge]
    y2_subgraph = f.edge_to_subgraph[g_nodes[2].i[:out].edge]
    y3_subgraph = f.edge_to_subgraph[g_nodes[3].i[:out].edge]

    @fact length(qs) => 9
    @fact qs[(g_nodes[1], m_subgraph)].distribution => vague(GaussianDistribution)
    @fact qs[(g_nodes[2], m_subgraph)].distribution => vague(GaussianDistribution)
    @fact qs[(g_nodes[3], m_subgraph)].distribution => vague(GaussianDistribution)
    @fact qs[(g_nodes[1], gam_subgraph)].distribution => vague(GammaDistribution)
    @fact qs[(g_nodes[2], gam_subgraph)].distribution => vague(GammaDistribution)
    @fact qs[(g_nodes[3], gam_subgraph)].distribution => vague(GammaDistribution)
    @fact qs[(g_nodes[1], y1_subgraph)].distribution => vague(GaussianDistribution)
    @fact qs[(g_nodes[2], y2_subgraph)].distribution => vague(GaussianDistribution)
    @fact qs[(g_nodes[3], y3_subgraph)].distribution => vague(GaussianDistribution)

    # Structured case
    data = [1.0]
    (g_nodes, y_nodes, m_eq_nodes, gam_eq_nodes, q_m_edges, q_gam_edges, q_y_edges) = initializeGaussianNodeChain(data)
    n_sections = length(data)

    f = VMP.QFactorization()
    for edge in q_y_edges
        f = VMP.factorize!(Set{Edge}({edge}), f)
    end
    VMP.generateSchedule!(f) # Generate and store internal and external schedules on factorization subgraphs
    qs = VMP.vagueQDistributions(f)

    m_gam_subgraph = f.edge_to_subgraph[g_nodes[1].i[:mean].edge]
    y1_subgraph = f.edge_to_subgraph[g_nodes[1].i[:out].edge]

    @fact length(qs) => 2
    @fact qs[(g_nodes[1], m_gam_subgraph)].distribution => vague(NormalGammaDistribution)
    @fact qs[(g_nodes[1], y1_subgraph)].distribution => vague(GaussianDistribution)
end