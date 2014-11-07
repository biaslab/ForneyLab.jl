#####################
# Integration tests
#####################

facts("Factorization integration tests") do
    context("factorize!()") do
        context("Should include argument edges in a new subgraph") do
            (driver, inhibitor, noise, add) = initializeLoopyGraph()
            factorize!(Set{Edge}({inhibitor.out.edge})) # Put this edge in a different subgraph
            graph = getCurrentGraph()
            @fact graph.factorization[1].nodes => Set{Node}({driver, inhibitor, noise, add})
            @fact graph.factorization[1].internal_edges => Set{Edge}({add.out.edge, add.in1.edge, add.in2.edge})
            @fact graph.factorization[1].external_edges => Set{Edge}({inhibitor.out.edge})
            @fact graph.factorization[2].nodes => Set{Node}({driver, inhibitor})
            @fact graph.factorization[2].internal_edges => Set{Edge}({inhibitor.out.edge})
            @fact graph.factorization[2].external_edges => Set{Edge}({add.in1.edge, add.out.edge})

            @fact true => false
        end

        context("Should update the edge_to_subgraph mapping for the graph") do
            (driver, inhibitor, noise, add) = initializeLoopyGraph()
            factorize!(Set{Edge}({inhibitor.out.edge})) # Put this edge in a different subgraph
            graph = getCurrentGraph()
            @fact graph.edge_to_subgraph[add.out.edge] => graph.factorization[1]
            @fact graph.edge_to_subgraph[add.in1.edge] => graph.factorization[1]
            @fact graph.edge_to_subgraph[add.in2.edge] => graph.factorization[1]
            @fact graph.edge_to_subgraph[inhibitor.out.edge] => graph.factorization[2]

            @fact true => false
        end

        context("Should extend the edge set to envelope deterministic nodes") do
            @fact true => false
        end
    end

    context("factorizeMeanField!() should output a mean field factorized graph") do
        data = [1.0, 1.0, 1.0]
        (g_nodes, y_nodes, m_eq_nodes, gam_eq_nodes, q_m_edges, q_gam_edges, q_y_edges) = initializeGaussianNodeChain(data)
        graph = getCurrentGraph()
        factorizeMeanField!(graph)
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
        @fact getSubgraph(g_nodes[1].mean.edge).internal_edges => m_set 
        @fact getSubgraph(g_nodes[1].precision.edge).internal_edges => gam_set 
        @fact getSubgraph(g_nodes[1].out.edge).internal_edges => Set{Edge}([q_y_edges[1]]) 
        @fact getSubgraph(g_nodes[2].out.edge).internal_edges => Set{Edge}([q_y_edges[2]])
        @fact getSubgraph(g_nodes[3].out.edge).internal_edges => Set{Edge}([q_y_edges[3]])
    end
end