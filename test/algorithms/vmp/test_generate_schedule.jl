facts("VMP.generateSchedule() tests") do
    context("Should include backward messages when there is only one internal interface connected to an external node") do
        data = [1.0]
        (g_node, y_node, m_0_node, gam_0_node, m_N_node, gam_N_node, m_eq_node, gam_eq_node, m_edge, gam_edge, y_edge) = initializeGaussianNodeChainForSvmp(data)

        # Structured factorization
        f = VMP.factorize!(Set{Edge}([y_edge]))

        graph = currentGraph()
        for subgraph in f.factors
            VMP.generateSchedule!(subgraph) # Generate internal and external schedule automatically
        end

        y_subgraph = f.edge_to_subgraph[y_edge]
        m_gam_subgraph = f.edge_to_subgraph[m_edge]
        @fact y_subgraph.internal_schedule => ForneyLab.convert(Schedule, [y_edge.head, y_edge.tail]) # Include outgoing interface
        @fact m_gam_subgraph.internal_schedule => ForneyLab.convert(Schedule, [m_0_node.out, m_N_node.out, m_edge.tail, gam_0_node.out, gam_N_node.out, gam_edge.tail]) # Exclude outgoing interfaces
    end
    context("Should generate an internal and external schedule when called on a subgraph") do
        (t1, a1, g1, t2, t3) = initializeFactoringGraphWithoutLoop()
        f = VMP.factorize!(Set{Edge}([t2.out.edge])) # Put this edge in a different subgraph
        for subgraph in f.factors
            VMP.generateSchedule!(subgraph)
            @fact length(unique(subgraph.internal_schedule)) => length(subgraph.internal_schedule) # No duplicate entries in schedule
        end
        # There are multiple valid schedules because of different orderings. Validity or schedule order is not checked here.
        @fact f.factors[1].internal_schedule => ForneyLab.convert(Schedule, [t1.out, a1.out, t3.out])
        @fact f.factors[2].internal_schedule => ForneyLab.convert(Schedule, [t2.out, t2.out.partner])
        @fact f.factors[1].external_schedule => [g1]
        @fact f.factors[2].external_schedule => [g1]
    end
end