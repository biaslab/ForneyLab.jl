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
            # All (but just) required calculations should be in the schedule
            @fact inhibitor.out in schedule => true
            @fact driver.out    in schedule => true
            @fact inhibitor.in1 in schedule => true
            @fact driver.in1    in schedule => true
            @fact add.in2       in schedule => true
            @fact add.in1       in schedule => false
            @fact add.out       in schedule => false
            @fact noise.out     in schedule => false
            # Validate correct relative order in schedule
            @fact findfirst(schedule, inhibitor.out)    < findfirst(schedule, driver.out)   => true
            @fact findfirst(schedule, driver.out)       < findfirst(schedule, add.in2)      => true
            @fact findfirst(schedule, driver.in1)       < findfirst(schedule, inhibitor.in1)=> true
            @fact findfirst(schedule, inhibitor.in1)    < findfirst(schedule, add.in2)      => true
        end

        context("Should correctly complete a partial schedule") do
            # Generate a schedule that first passes clockwise through the cycle and then counterclockwise
            schedule = generateSchedule([driver.out, add.in2]) # Message towards noise factor
            # All (but just) required calculations should be in the schedule
            @fact schedule[1] => inhibitor.out
            @fact schedule[2] => driver.out
            @fact schedule[3] => driver.in1
            @fact schedule[4] => inhibitor.in1
            @fact schedule[5] => add.in2
        end


        context("Should generate an internal and external schedule when called on a subgraph") do
            (driver, inhibitor, noise, add) = initializeLoopyGraph()
            factorize!(Set{Edge}({inhibitor.out.edge})) # Put this edge in a different subgraph
            graph = getCurrentGraph()
            for subgraph in graph.factorization
                generateSchedule!(subgraph)
                @fact length(unique(subgraph.internal_schedule)) => length(subgraph.internal_schedule) # No duplicate entries in schedule
            end
            # There are multiple valid schedules because of different orderings. Validity or schedule order is not checked here.
            @fact Set{Interface}(graph.factorization[1].internal_schedule) => Set{Interface}([noise.out, inhibitor.in1, add.in1, driver.out, add.out])
            @fact Set{Interface}(graph.factorization[2].internal_schedule) => Set{Interface}([driver.in1, inhibitor.out])
            @fact Set{Node}(graph.factorization[1].external_schedule) => Set{Node}([driver, inhibitor])
            @fact Set{Node}(graph.factorization[2].external_schedule) => Set{Node}([inhibitor, driver])
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
            @fact Set(y_subgraph.internal_schedule) => Set([y_edge.head, y_edge.tail])#Include outgoing interface
            @fact Set(m_gam_subgraph.internal_schedule) => Set([m_edge.tail, gam_edge.tail, m_0_node.out, gam_0_node.out, m_N_node.out, gam_N_node.out])# Exclude outgoing interfaces
        end
    end
end