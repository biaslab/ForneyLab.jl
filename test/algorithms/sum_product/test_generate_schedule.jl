facts("Schedule related tests") do
    context("ForneyLab.convert(Schedule, ...)") do
        node = GaussianNode()
        @fact ForneyLab.convert(Schedule, [node.out, node.mean]) => [ScheduleEntry(node.out, sumProduct!), ScheduleEntry(node.mean, sumProduct!)]
        @fact ForneyLab.convert(Schedule, [node.out, node.mean], sumProduct!, sample) => [ScheduleEntry(node.out, sumProduct!, sample), ScheduleEntry(node.mean, sumProduct!, sample)]
    end
end

facts("SumProduct.generateSchedule() integration tests") do
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
            schedule = SumProduct.generateSchedule(add.in2) # Message towards noise factor
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
            schedule = SumProduct.generateSchedule(ForneyLab.convert(Schedule, [driver.out, add.in2])) # Message towards noise factor
            # All (but just) required calculations should be in the schedule
            @fact schedule => ForneyLab.convert(Schedule, [inhibitor.out, driver.out, driver.in1, inhibitor.in1, add.in2])
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
            @fact s2 => ForneyLab.convert(Schedule, [node_t1.out.partner])
        end

        context("Should generate a schedule that propagates messages to write buffers defined on interfaces") do
            g = FactorGraph()
            node_t1 = TerminalNode()
            node_t2 = TerminalNode()
            e = Edge(node_t1, node_t2)
            s1 = SumProduct.generateSchedule(g)
            @fact s1 => Array(ScheduleEntry, 0)
            setWriteBuffer(node_t1.out)
            s2 = SumProduct.generateSchedule(g)
            @fact s2 => ForneyLab.convert(Schedule, [node_t1.out])
        end

        context("Should generate a schedule that propagates messages to write buffers defined on edges") do
            g = FactorGraph()
            node_t1 = TerminalNode()
            node_t2 = TerminalNode()
            e = Edge(node_t1, node_t2)
            s1 = SumProduct.generateSchedule(g)
            @fact s1 => Array(ScheduleEntry, 0)
            setWriteBuffer(node_t1.out.edge)
            s2 = SumProduct.generateSchedule(g)
            @fact s2 => ForneyLab.convert(Schedule, [node_t1.out.partner, node_t1.out])
        end
    end
end
