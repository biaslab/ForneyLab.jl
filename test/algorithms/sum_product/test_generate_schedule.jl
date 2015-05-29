facts("Schedule related tests") do
    context("ForneyLab.convert(Schedule, ...)") do
        node = GaussianNode()
        @fact ForneyLab.convert(Schedule, [node.i[:out], node.i[:mean]]) => [ScheduleEntry(node.i[:out], sumProduct!), ScheduleEntry(node.i[:mean], sumProduct!)]
        @fact ForneyLab.convert(Schedule, [node.i[:out], node.i[:mean]], sumProduct!, sample) => [ScheduleEntry(node.i[:out], sumProduct!, sample), ScheduleEntry(node.i[:mean], sumProduct!, sample)]
    end
end

facts("SumProduct.generateSchedule() integration tests") do
    context("generateSchedule()") do
        (driver, inhibitor, noise, add) = initializeLoopyGraph(A=[2.0], B=[0.5], noise_m=1.0, noise_V=0.1)

        context("Should throw an error when there is an unbroken loop") do
            @fact_throws generateSchedule(driver.i[:out])
        end

        # Initial message
        setMessage!(add.i[:in1], Message(GaussianDistribution(m=2.0, V=0.5)))
        setMessage!(add.i[:out], Message(GaussianDistribution()))

        context("Should auto-generate a feasible schedule") do
            # Generate schedule automatically
            schedule = SumProduct.generateSchedule(add.i[:in2]) # Message towards noise factor
            intf_list = [schedule_entry.interface for schedule_entry in schedule]
            # All (but just) required calculations should be in the schedule
            @fact inhibitor.i[:out] in intf_list => true
            @fact driver.i[:out]    in intf_list => true
            @fact inhibitor.i[:in] in intf_list => true
            @fact driver.i[:in]    in intf_list => true
            @fact add.i[:in2]       in intf_list => true
            @fact add.i[:in1]       in intf_list => false
            @fact add.i[:out]       in intf_list => false
            @fact noise.i[:out]     in intf_list => false
            # Validate correct relative order in schedule
            @fact findfirst(intf_list, inhibitor.i[:out])    < findfirst(intf_list, driver.i[:out])   => true
            @fact findfirst(intf_list, driver.i[:out])       < findfirst(intf_list, add.i[:in2])      => true
            @fact findfirst(intf_list, driver.i[:in])       < findfirst(intf_list, inhibitor.i[:in])=> true
            @fact findfirst(intf_list, inhibitor.i[:in])    < findfirst(intf_list, add.i[:in2])      => true
        end

        context("Should correctly complete a partial schedule") do
            # Generate a schedule that first passes clockwise through the cycle and then counterclockwise
            schedule = SumProduct.generateSchedule(ForneyLab.convert(Schedule, [driver.i[:out], add.i[:in2]])) # Message towards noise factor
            # All (but just) required calculations should be in the schedule
            @fact schedule => ForneyLab.convert(Schedule, [inhibitor.i[:out], driver.i[:out], driver.i[:in], inhibitor.i[:in], add.i[:in2]])
        end

        context("Should generate a schedule that propagates messages to timewraps") do
            g = FactorGraph()
            node_t1 = TerminalNode()
            node_t2 = TerminalNode()
            e = Edge(node_t1, node_t2)
            s1 = SumProduct.generateSchedule(g)
            @fact s1 => Array(ScheduleEntry, 0)
            wrap(node_t1, node_t2)
            s2 = SumProduct.generateSchedule(g)
            @fact s2 => ForneyLab.convert(Schedule, [node_t1.i[:out].partner])
        end

        context("Should generate a schedule that propagates messages to write buffers defined on interfaces") do
            g = FactorGraph()
            node_t1 = TerminalNode()
            node_t2 = TerminalNode()
            e = Edge(node_t1, node_t2)
            s1 = SumProduct.generateSchedule(g)
            @fact s1 => Array(ScheduleEntry, 0)
            setWriteBuffer(node_t1.i[:out])
            s2 = SumProduct.generateSchedule(g)
            @fact s2 => ForneyLab.convert(Schedule, [node_t1.i[:out]])
        end

        context("Should generate a schedule that propagates messages to write buffers defined on edges") do
            g = FactorGraph()
            node_t1 = TerminalNode()
            node_t2 = TerminalNode()
            e = Edge(node_t1, node_t2)
            s1 = SumProduct.generateSchedule(g)
            @fact s1 => Array(ScheduleEntry, 0)
            setWriteBuffer(node_t1.i[:out].edge)
            s2 = SumProduct.generateSchedule(g)
            @fact s2 => ForneyLab.convert(Schedule, [node_t1.i[:out].partner, node_t1.i[:out]])
        end
    end
end
