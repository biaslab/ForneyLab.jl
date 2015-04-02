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
            scheme = InferenceScheme()
            factorize!(Set{Edge}([t2.out.edge]))
            graph = currentGraph()
            @fact scheme.factorization[2].internal_edges => Set{Edge}({t2.out.edge})
            @fact scheme.factorization[1].internal_edges => Set{Edge}({t1.out.edge, a1.out.edge, g1.out.edge, add1.out.edge, g2.out.edge})
        end

        context("Should update the edge_to_subgraph mapping for the graph") do
            (t1, a1, g1, t2, add1, g2) = initializeFactoringGraph()
            scheme = InferenceScheme()
            factorize!(Set{Edge}([t2.out.edge]))
            graph = currentGraph()
            @fact scheme.edge_to_subgraph[t1.out.edge] => scheme.factorization[1]
            @fact scheme.edge_to_subgraph[a1.out.edge] => scheme.factorization[1]
            @fact scheme.edge_to_subgraph[g1.out.edge] => scheme.factorization[1]
            @fact scheme.edge_to_subgraph[t2.out.edge] => scheme.factorization[2]
            @fact scheme.edge_to_subgraph[add1.out.edge] => scheme.factorization[1]
            @fact scheme.edge_to_subgraph[g2.out.edge] => scheme.factorization[1]
        end

        context("Should not factorize a GaussianNode with fixed mean and variance") do
            (t, gain, gauss) = initializeGaussianFactoringGraph()
            scheme = InferenceScheme()
            factorize!(Set{Edge}([t.out.edge]))
            graph = currentGraph()
            @fact length(scheme.factorization) => 1
            @fact scheme.edge_to_subgraph[t.out.edge] => scheme.factorization[1]
            @fact scheme.edge_to_subgraph[gauss.out.edge] => scheme.factorization[1]
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
        scheme = currentScheme()
        factorize!(scheme)
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
        @fact length(scheme.factorization) => 5

        @fact subgraph(scheme, g_nodes[1].mean.edge).internal_edges => m_set 
        @fact subgraph(scheme, g_nodes[1].precision.edge).internal_edges => gam_set 
        @fact subgraph(scheme, g_nodes[1].out.edge).internal_edges => Set{Edge}([q_y_edges[1]]) 
        @fact subgraph(scheme, g_nodes[2].out.edge).internal_edges => Set{Edge}([q_y_edges[2]])
        @fact subgraph(scheme, g_nodes[3].out.edge).internal_edges => Set{Edge}([q_y_edges[3]])
    end
end