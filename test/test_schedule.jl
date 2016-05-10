facts("Schedule and ScheduleEntry tests") do
    context("General properties") do
        FactorGraph()

        # ScheduleEntry
        entry1 = ScheduleEntry(MockNode(id=:mock1), 1, sumProductRule!)
        entry2 = ScheduleEntry(MockNode(id=:mock2), 1, sumProductRule!)
        @fact is(entry1.node.interfaces[entry1.outbound_interface_id], n(:mock1).i[:out]) --> true
        @fact_throws deepcopy(entry1)

        # copy(::ScheduleEntry)
        entry1_copy = copy(entry1)
        @fact is(entry1, entry1_copy) --> false
        @fact is(entry1.node, entry1_copy.node) --> true
        @fact is(entry1.rule, entry1_copy.rule) --> true
        @fact isdefined(entry1_copy, :involves_approximation) --> isdefined(entry1, :involves_approximation)
        if isdefined(entry1_copy, :involves_approximation)
            @fact entry1_copy.involves_approximation --> entry1.involves_approximation
        end

        # Schedule
        schedule = [entry1, entry2]
        @fact typeof(schedule) --> Schedule
        schedule_deepcopy = deepcopy(schedule)
        for idx=1:length(schedule)
            @fact is(schedule[idx], schedule_deepcopy[idx]) --> false
            @fact is(schedule[idx].node, schedule_deepcopy[idx].node) --> true
        end
    end

    context("A compiled scheduleEntry can be executed") do
        node = TerminalNode(Gaussian(m=2.0, V=4.0))
        node.i[:out].message = Message(Gaussian()) # Preset message
        entry = ScheduleEntry(node, 1, sumProductRule!)
        ForneyLab.buildExecute!(entry, Any[nothing])
        outbound_dist = entry.execute()
        @fact is(node.i[:out].message.payload, outbound_dist) --> true
        @fact node.i[:out].message --> Message(Gaussian(m=2.0, V=4.0))

        node = AdditionNode()
        node.i[:out].message = Message(Gaussian()) # Preset message
        entry = ScheduleEntry(node, 3, sumProductRule!)
        ForneyLab.buildExecute!(entry, Any[Message(Gaussian(m=1.0, V=1.0)), Message(Gaussian(m=1.0, V=1.0)), nothing])
        outbound_dist = entry.execute()
        @fact is(node.i[:out].message.payload, outbound_dist) --> true
        @fact node.i[:out].message --> Message(Gaussian(m=2.0, V=2.0))

        node = GaussianNode(form=:precision)
        node.i[:out].message = Message(Gaussian()) # Preset message
        entry = ScheduleEntry(node, 3, variationalRule!)
        ForneyLab.buildExecute!(entry, Any[Gaussian(m=1.0, V=1.0), Gamma(a=2.0, b=4.0), nothing])
        outbound_dist = entry.execute()
        @fact is(node.i[:out].message.payload, outbound_dist) --> true
        @fact node.i[:out].message --> Message(Gaussian(m=1.0, W=0.5))
    end
end
