facts("Schedule related tests") do
    context("ForneyLab.convert(Schedule, ...)") do
        FactorGraph()
        GaussianNode(id=:node)
        @fact ForneyLab.convert(Schedule, [n(:node).i[:out], n(:node).i[:mean]]) => [ScheduleEntry(n(:node).i[:out], sumProduct!), ScheduleEntry(n(:node).i[:mean], sumProduct!)]
        @fact ForneyLab.convert(Schedule, [n(:node).i[:out], n(:node).i[:mean]], sumProduct!, sample) => [ScheduleEntry(n(:node).i[:out], sumProduct!, sample), ScheduleEntry(n(:node).i[:mean], sumProduct!, sample)]
    end
end

facts("SumProduct.generateSchedule() integration tests") do
    context("generateSchedule()") do
        initializeLoopyGraph(A=[2.0], B=[0.5], noise_m=1.0, noise_V=0.1)

        context("Should throw an error when there is an unbroken loop") do
            @fact_throws generateSchedule(n(:driver).i[:out])
        end

        # Initial message
        setMessage!(n(:add).i[:in1], Message(GaussianDistribution(m=2.0, V=0.5)))
        setMessage!(n(:add).i[:out], Message(GaussianDistribution()))

        context("Should auto-generate a feasible schedule") do
            # Generate schedule automatically
            schedule = SumProduct.generateSchedule(n(:add).i[:in2]) # Message towards noise factor
            intf_list = [schedule_entry.interface for schedule_entry in schedule]
            # All (but just) required calculations should be in the schedule
            @fact n(:inhibitor).i[:out] in intf_list => true
            @fact n(:driver).i[:out]    in intf_list => true
            @fact n(:inhibitor).i[:in] in intf_list => true
            @fact n(:driver).i[:in]    in intf_list => true
            @fact n(:add).i[:in2]       in intf_list => true
            @fact n(:add).i[:in1]       in intf_list => false
            @fact n(:add).i[:out]       in intf_list => false
            @fact n(:noise).i[:out]     in intf_list => false
            # Validate correct relative order in schedule
            @fact findfirst(intf_list, n(:inhibitor).i[:out])    < findfirst(intf_list, n(:driver).i[:out])   => true
            @fact findfirst(intf_list, n(:driver).i[:out])       < findfirst(intf_list, n(:add).i[:in2])      => true
            @fact findfirst(intf_list, n(:driver).i[:in])       < findfirst(intf_list, n(:inhibitor).i[:in])=> true
            @fact findfirst(intf_list, n(:inhibitor).i[:in])    < findfirst(intf_list, n(:add).i[:in2])      => true
        end

        context("Should correctly complete a partial schedule") do
            # Generate a schedule that first passes clockwise through the cycle and then counterclockwise
            schedule = SumProduct.generateSchedule(ForneyLab.convert(Schedule, [n(:driver).i[:out], n(:add).i[:in2]])) # Message towards noise factor
            # All (but just) required calculations should be in the schedule
            @fact schedule => ForneyLab.convert(Schedule, [n(:inhibitor).i[:out], n(:driver).i[:out], n(:driver).i[:in], n(:inhibitor).i[:in], n(:add).i[:in2]])
        end

        context("Should generate a schedule that propagates messages to timewraps") do
            g = FactorGraph()
            TerminalNode(id=:t1)
            TerminalNode(id=:t2)
            e = Edge(n(:t1), n(:t2))
            s1 = SumProduct.generateSchedule(g)
            @fact s1 => Array(ScheduleEntry, 0)
            wrap(n(:t1), n(:t2))
            s2 = SumProduct.generateSchedule(g)
            @fact s2 => ForneyLab.convert(Schedule, [n(:t1).i[:out].partner])
        end

        context("Should generate a schedule that propagates messages to write buffers defined on interfaces") do
            g = FactorGraph()
            TerminalNode(id=:t1)
            TerminalNode(id=:t2)
            e = Edge(n(:t1), n(:t2))
            s1 = SumProduct.generateSchedule(g)
            @fact s1 => Array(ScheduleEntry, 0)
            setWriteBuffer(n(:t1).i[:out])
            s2 = SumProduct.generateSchedule(g)
            @fact s2 => ForneyLab.convert(Schedule, [n(:t1).i[:out]])
        end

        context("Should generate a schedule that propagates messages to write buffers defined on edges") do
            g = FactorGraph()
            TerminalNode(id=:t1)
            TerminalNode(id=:t2)
            e = Edge(n(:t1), n(:t2))
            s1 = SumProduct.generateSchedule(g)
            @fact s1 => Array(ScheduleEntry, 0)
            setWriteBuffer(n(:t1).i[:out].edge)
            s2 = SumProduct.generateSchedule(g)
            @fact s2 => ForneyLab.convert(Schedule, [n(:t1).i[:out].partner, n(:t1).i[:out]])
        end
    end
end
