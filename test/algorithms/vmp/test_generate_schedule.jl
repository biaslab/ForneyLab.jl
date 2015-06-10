facts("VMP.generateSchedule() tests") do
    context("Should include backward messages when there is only one internal interface connected to an external node and include vmp updates") do
        data = [1.0]
        initializeGaussianNodeChainForSvmp(data)

        # Structured factorization
        f = VMP.factorize!(Set{Edge}([e(:y)]))

        graph = currentGraph()
        for subgraph in f.factors
            VMP.generateSchedule!(subgraph) # Generate internal and external schedule automatically
        end

        y_subgraph = f.edge_to_subgraph[e(:y)]
        m_gam_subgraph = f.edge_to_subgraph[e(:m)]
        @fact y_subgraph.internal_schedule => [ForneyLab.ScheduleEntry(e(:y).head, sumProduct!), ForneyLab.ScheduleEntry(e(:y).tail, vmp!)] # Include outgoing interface
        @fact m_gam_subgraph.internal_schedule => ForneyLab.convert(Schedule, [n(:m0).i[:out], n(:mN).i[:out], e(:m).tail, n(:gam0).i[:out], n(:gamN).i[:out], e(:gam).tail]) # Exclude outgoing interfaces
    end
    context("Should generate an internal and external schedule when called on a subgraph") do
        initializeFactoringGraphWithoutLoop()
        f = VMP.factorize!(Set{Edge}([n(:t2).i[:out].edge])) # Put this edge in a different subgraph
        for subgraph in f.factors
            VMP.generateSchedule!(subgraph)
            @fact length(unique(subgraph.internal_schedule)) => length(subgraph.internal_schedule) # No duplicate entries in schedule
        end
        # There are multiple valid schedules because of different orderings. Validity or schedule order is not checked here.
        @fact f.factors[1].internal_schedule => ForneyLab.convert(Schedule, [n(:t1).i[:out], n(:a1).i[:out], n(:t3).i[:out]])
        @fact f.factors[2].internal_schedule => [ForneyLab.ScheduleEntry(n(:t2).i[:out], sumProduct!), ForneyLab.ScheduleEntry(n(:t2).i[:out].partner, vmp!)]
        @fact f.factors[1].external_schedule => [n(:g1)]
        @fact f.factors[2].external_schedule => [n(:g1)]
    end
end