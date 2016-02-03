facts("ForneyLab.generateVariationalBayesSchedule() tests") do
    context("Should include backward messages when there is only one internal interface connected to an external node and include vmp updates") do
        data = [1.0]
        initializeGaussianNodeChain(data)

        # Structured factorization
        f = ForneyLab.factorize!(Set{Edge}([edge(:q_y1)]))

        graph = currentGraph()
        for subgraph in f.factors
            ForneyLab.generateVariationalBayesSchedule!(subgraph) # Generate internal and external schedule automatically
        end

        y_subgraph = f.edge_to_subgraph[edge(:q_y1)]
        m_gam_subgraph = f.edge_to_subgraph[edge(:q_m1)]
        @fact y_subgraph.internal_schedule --> [ForneyLab.ScheduleEntry(edge(:q_y1).head, sumProductRule!), ForneyLab.ScheduleEntry(edge(:q_y1).tail, variationalRule!)] # Include outgoing interface
        @fact m_gam_subgraph.internal_schedule --> ForneyLab.convert(Schedule, [n(:m0).i[:out], n(:mN).i[:out], edge(:q_m1).tail, n(:gam0).i[:out], n(:gamN).i[:out], eg(:q_gam1).tail]) # Exclude outgoing interfaces
    end
    context("Should generate an internal and external schedule when called on a subgraph") do
        initializeFactoringGraphWithoutLoop()
        f = ForneyLab.factorize!(Set{Edge}([n(:t2).i[:out].edge])) # Put this edge in a different subgraph
        for subgraph in f.factors
            ForneyLab.generateVariationalBayesSchedule!(subgraph)
            @fact length(unique(subgraph.internal_schedule)) --> length(subgraph.internal_schedule) # No duplicate entries in schedule
        end
        # There are multiple valid schedules because of different orderings. Validity or schedule order is not checked here.
        @fact f.factors[1].internal_schedule --> ForneyLab.convert(Schedule, [n(:t1).i[:out], n(:a1).i[:out], n(:t3).i[:out]])
        @fact f.factors[2].internal_schedule --> [ForneyLab.ScheduleEntry(n(:t2).i[:out], sumProductRule!), ForneyLab.ScheduleEntry(n(:t2).i[:out].partner, variationalRule!)]
        @fact f.factors[1].external_schedule --> [n(:g1)]
        @fact f.factors[2].external_schedule --> [n(:g1)]
    end
end
