facts("Schedule related tests") do
    context("ForneyLab.convert(Schedule, ...)") do
        node = GaussianNode()
        @fact ForneyLab.convert(Schedule, [node.out, node.mean]) => [ScheduleEntry(node.out, :sumproduct), ScheduleEntry(node.mean, :sumproduct)]
        @fact ForneyLab.convert(Schedule, [node.out, node.mean], :sumproduct_sample) => [ScheduleEntry(node.out, :sumproduct_sample), ScheduleEntry(node.mean, :sumproduct_sample)]
        @fact ForneyLab.convert(Schedule, [node.out, node.mean], :sumproduct_expectation) => [ScheduleEntry(node.out, :sumproduct_expectation), ScheduleEntry(node.mean, :sumproduct_expectation)]
    end
end

facts("generateSchedule() integration tests") do
    context("generateSchedule()") do
        (driver, inhibitor, noise, add) = initializeLoopyGraph(A=[2.0], B=[0.5], noise_m=1.0, noise_V=0.1)

        context("Should throw an error when there is an unbroken loop") do
            @fact_throws generateSchedule(driver.out)
        end

        # Initial message
        setMessage!(add.in1, Message(GaussianDistribution(m=2.0, V=0.5)))
        setMessage!(add.out, Message(GaussianDistribution()))

        context("Should auto-generate a feasible schedule") do
            # Generate schedule automatically
            schedule = generateSchedule(add.in2) # Message towards noise factor
            intf_list = [schedule_entry.interface for schedule_entry in schedule]
            # All (but just) required calculations should be in the schedule
            @fact inhibitor.out in intf_list => true
            @fact driver.out    in intf_list => true
            @fact inhibitor.in1 in intf_list => true
            @fact driver.in1    in intf_list => true
            @fact add.in2       in intf_list => true
            @fact add.in1       in intf_list => false
            @fact add.out       in intf_list => false
            @fact noise.out     in intf_list => false
            # Validate correct relative order in schedule
            @fact findfirst(intf_list, inhibitor.out)    < findfirst(intf_list, driver.out)   => true
            @fact findfirst(intf_list, driver.out)       < findfirst(intf_list, add.in2)      => true
            @fact findfirst(intf_list, driver.in1)       < findfirst(intf_list, inhibitor.in1)=> true
            @fact findfirst(intf_list, inhibitor.in1)    < findfirst(intf_list, add.in2)      => true
        end

        context("Should correctly complete a partial schedule") do
            # Generate a schedule that first passes clockwise through the cycle and then counterclockwise
            schedule = generateSchedule(ForneyLab.convert(Schedule, [driver.out, add.in2])) # Message towards noise factor
            # All (but just) required calculations should be in the schedule
            @fact schedule => ForneyLab.convert(Schedule, [inhibitor.out, driver.out, driver.in1, inhibitor.in1, add.in2])
        end

        context("Should generate an internal and external schedule when called on a subgraph") do
            (t1, a1, g1, t2, t3) = initializeFactoringGraphWithoutLoop()
            factorize!(Set{Edge}([t2.out.edge])) # Put this edge in a different subgraph
            graph = getCurrentGraph()
            for subgraph in graph.factorization
                generateSchedule!(subgraph)
                @fact length(unique(subgraph.internal_schedule)) => length(subgraph.internal_schedule) # No duplicate entries in schedule
            end
            # There are multiple valid schedules because of different orderings. Validity or schedule order is not checked here.
            @fact graph.factorization[1].internal_schedule => ForneyLab.convert(Schedule, [t1.out, a1.out, t3.out])
            @fact graph.factorization[2].internal_schedule => ForneyLab.convert(Schedule, [t2.out, t2.out.partner])
            @fact graph.factorization[1].external_schedule => [g1]
            @fact graph.factorization[2].external_schedule => [g1]
        end

        context("Should generate a schedule that propagates messages to timewraps when called on a subgraph") do
            g = FactorGraph()
            node_t1 = TerminalNode()
            node_t2 = TerminalNode()
            e = Edge(node_t1, node_t2)
            generateSchedule!(g.factorization[1])
            @fact g.factorization[1].internal_schedule => Array(ScheduleEntry, 0)
            addTimeWrap(node_t1, node_t2)
            generateSchedule!(g.factorization[1])
            @fact g.factorization[1].internal_schedule => ForneyLab.convert(Schedule, [node_t1.out.partner])
        end

        context("Should include backward messages when there is only one internal interface connected to an external node") do
            data = [1.0]
            (g_node, y_node, m_0_node, gam_0_node, m_N_node, gam_N_node, m_eq_node, gam_eq_node, m_edge, gam_edge, y_edge) = initializeGaussianNodeChainForSvmp(data)

            # Structured factorization
            factorize!(Set{Edge}([y_edge]))

            graph = getCurrentGraph()
            for subgraph in graph.factorization
                generateSchedule!(subgraph) # Generate internal and external schedule automatically
            end

            y_subgraph = getSubgraph(y_edge)
            m_gam_subgraph = getSubgraph(m_edge)
            @fact y_subgraph.internal_schedule => ForneyLab.convert(Schedule, [y_edge.head, y_edge.tail]) # Include outgoing interface
            @fact m_gam_subgraph.internal_schedule => ForneyLab.convert(Schedule, [m_0_node.out, m_N_node.out, m_edge.tail, gam_0_node.out, gam_N_node.out, gam_edge.tail]) # Exclude outgoing interfaces
        end
    end
end