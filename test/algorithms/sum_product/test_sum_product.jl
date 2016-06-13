facts("SumProduct") do
    context("SumProduct should perform message type inference") do
        initializeAdditionNode([Gaussian(), Gaussian(), Gaussian()])

        # Forcing a message to a distribution type for which no rule is implemented should throw an error
        message_types = Dict{Interface,DataType}(n(:add_node).i[:out] => Gamma)
        @fact_throws SumProduct(n(:add_node).i[:out], message_types=message_types)

        # Forcing a message to a supported distribution type should be ok.
        message_types = Dict{Interface,DataType}(n(:add_node).i[:out] => Gaussian)
        algo = SumProduct(n(:add_node).i[:out], message_types=message_types)
        @fact algo.schedule[end].inbound_types --> [Message{Gaussian}, Message{Gaussian}, Void]
        @fact algo.schedule[end].outbound_type --> Gaussian
        @fact isdefined(algo.schedule[end], :approximation) --> false
    end

    context("Should generate a schedule that propagates messages to wraps") do
        g = FactorGraph()
        TerminalNode(id=:t1)
        TerminalNode(id=:t2)
        e = Edge(n(:t1), n(:t2))
        algo = SumProduct(g)
        @fact length(algo.schedule) --> 0
        Wrap(n(:t1), n(:t2))
        algo2 = SumProduct(g)
        @fact algo2.schedule --> ForneyLab.convert(Schedule, [n(:t1).i[:out].partner], SumProductRule)
    end

    context("Should generate a schedule that propagates messages to write buffers defined on interfaces") do
        g = FactorGraph()
        TerminalNode(id=:t1)
        TerminalNode(id=:t2)
        e = Edge(n(:t1), n(:t2))
        algo = SumProduct(g)
        @fact length(algo.schedule) --> 0
        attachWriteBuffer(n(:t1).i[:out])
        algo = SumProduct(g)
        @fact algo.schedule --> ForneyLab.convert(Schedule, [n(:t1).i[:out]], SumProductRule)
    end

    context("Should generate a schedule that propagates messages to write buffers defined on edges") do
        g = FactorGraph()
        TerminalNode(id=:t1)
        TerminalNode(id=:t2)
        e = Edge(n(:t1), n(:t2))
        algo = SumProduct(g)
        @fact length(algo.schedule) --> 0
        attachWriteBuffer(n(:t1).i[:out].edge)
        algo = SumProduct(g)
        @fact algo.schedule --> ForneyLab.convert(Schedule, [n(:t1).i[:out].partner, n(:t1).i[:out]], SumProductRule)
    end
end


facts("SumProduct message passing tests") do
    context("SumProduct execute()") do
        context("Should correctly execute a schedule and return the result of the last step") do
            initializeAdditionNode()

            algo = SumProduct(n(:add_node).i[:out])
            prepare!(algo)
            msg = execute(algo)

            @fact msg.payload --> Gaussian(m=0.0, V=2.0)
        end
    end
end
