#####################
# Integration tests
#####################

facts("Factorization integration tests") do
    context("extend() should extend a set of edges to envelope deterministic nodes") do
        (t1, a1, g1, t2, add1, g2) = initializeFactoringGraph()
        cluster = ForneyLab.extend(Set{Edge}([g1.out.edge, g1.mean.edge]))
        @fact cluster => Set{Edge}({t1.out.edge, a1.out.edge, g1.out.edge, add1.out.edge, g2.out.edge})
    end

    context("factorize!()") do
        context("Should include argument edges in a new subgraph") do
            (t1, a1, g1, t2, add1, g2) = initializeFactoringGraph()
            factorize!(Set{Edge}([t2.out.edge]))
            graph = currentGraph()
            @fact graph.factorization[2].nodes => Set{Node}({t2, g1})
            @fact graph.factorization[2].internal_edges => Set{Edge}({t2.out.edge})
            @fact graph.factorization[2].external_edges => Set{Edge}({g1.mean.edge, g1.out.edge})
            @fact graph.factorization[1].nodes => Set{Node}({t1, a1, g1, add1, g2})
            @fact graph.factorization[1].internal_edges => Set{Edge}({t1.out.edge, a1.out.edge, g1.out.edge, add1.out.edge, g2.out.edge})
            @fact graph.factorization[1].external_edges => Set{Edge}({t2.out.edge})
        end

        context("Should update the edge_to_subgraph mapping for the graph") do
            (t1, a1, g1, t2, add1, g2) = initializeFactoringGraph()
            factorize!(Set{Edge}([t2.out.edge]))
            graph = currentGraph()
            @fact graph.edge_to_subgraph[t1.out.edge] => graph.factorization[1]
            @fact graph.edge_to_subgraph[a1.out.edge] => graph.factorization[1]
            @fact graph.edge_to_subgraph[g1.out.edge] => graph.factorization[1]
            @fact graph.edge_to_subgraph[t2.out.edge] => graph.factorization[2]
            @fact graph.edge_to_subgraph[add1.out.edge] => graph.factorization[1]
            @fact graph.edge_to_subgraph[g2.out.edge] => graph.factorization[1]
        end

        context("Should not factorize a GaussianNode with fixed mean and variance") do
            (t, gain, gauss) = initializeGaussianFactoringGraph()
            factorize!(Set{Edge}([t.out.edge]))
            graph = currentGraph()
            @fact length(graph.factorization) => 1
            @fact graph.edge_to_subgraph[t.out.edge] => graph.factorization[1]
            @fact graph.edge_to_subgraph[gauss.out.edge] => graph.factorization[1]
        end

        context("Should retain node names") do
            (t1, gain, t2) = initializeSimpleFactoringGraph()
            factorize!(Set{Edge}([t1.out.edge]))
            graph = currentGraph()
            @fact node("t1") => t1
        end
    end

    context("factorize!() should output a mean field factorized graph") do
        data = [1.0, 1.0, 1.0]
        (g_nodes, y_nodes, m_eq_nodes, gam_eq_nodes, q_m_edges, q_gam_edges, q_y_edges) = initializeGaussianNodeChain(data)
        graph = currentGraph()
        factorize!(graph)
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
        @fact length(graph.factorization) => 5

        @fact subgraph(g_nodes[1].mean.edge).internal_edges => m_set 
        @fact subgraph(g_nodes[1].precision.edge).internal_edges => gam_set 
        @fact subgraph(g_nodes[1].out.edge).internal_edges => Set{Edge}([q_y_edges[1]]) 
        @fact subgraph(g_nodes[2].out.edge).internal_edges => Set{Edge}([q_y_edges[2]])
        @fact subgraph(g_nodes[3].out.edge).internal_edges => Set{Edge}([q_y_edges[3]])

        @fact subgraph(g_nodes[1].mean.edge).external_edges => Set{Edge}([g_nodes[1].out.edge, g_nodes[1].precision.edge, g_nodes[2].out.edge, g_nodes[2].precision.edge, g_nodes[3].out.edge, g_nodes[3].precision.edge])
        @fact subgraph(g_nodes[1].precision.edge).external_edges => Set{Edge}([g_nodes[1].out.edge, g_nodes[1].mean.edge, g_nodes[2].out.edge, g_nodes[2].mean.edge, g_nodes[3].out.edge, g_nodes[3].mean.edge])
        @fact subgraph(g_nodes[1].out.edge).external_edges => Set{Edge}([g_nodes[1].mean.edge, g_nodes[1].precision.edge])
        @fact subgraph(g_nodes[2].out.edge).external_edges => Set{Edge}([g_nodes[2].mean.edge, g_nodes[2].precision.edge])
        @fact subgraph(g_nodes[3].out.edge).external_edges => Set{Edge}([g_nodes[3].mean.edge, g_nodes[3].precision.edge])

        @fact subgraph(g_nodes[1].mean.edge).nodes => Set{Node}([m_eq_nodes, g_nodes, m_eq_nodes[1].interfaces[1].partner.node, m_eq_nodes[end].interfaces[2].partner.node])
        @fact subgraph(g_nodes[1].precision.edge).nodes => Set{Node}([gam_eq_nodes, g_nodes, gam_eq_nodes[1].interfaces[1].partner.node, gam_eq_nodes[end].interfaces[2].partner.node])
        @fact subgraph(g_nodes[1].out.edge).nodes => Set{Node}([g_nodes[1], g_nodes[1].out.partner.node]) 
        @fact subgraph(g_nodes[2].out.edge).nodes => Set{Node}([g_nodes[2], g_nodes[2].out.partner.node])
        @fact subgraph(g_nodes[3].out.edge).nodes => Set{Node}([g_nodes[3], g_nodes[3].out.partner.node])
    end
end