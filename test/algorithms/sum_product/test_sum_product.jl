facts("SumProduct collect inbound type tests") do
    context("SumProduct should perform message type inference") do
        initializeAdditionNode([GaussianDistribution(), GaussianDistribution(), GaussianDistribution()])

        # Forcing a message to a distribution type for which no rule is implemented should throw an error
        message_types = Dict{Interface,DataType}(n(:add_node).i[:out] => GammaDistribution)
        @fact_throws SumProduct(n(:add_node).i[:out], message_types=message_types)

        # Forcing a message to a supported distribution type should be ok.
        message_types = Dict{Interface,DataType}(n(:add_node).i[:out] => GaussianDistribution)
        algo = SumProduct(n(:add_node).i[:out], message_types=message_types)
        @fact algo.schedule[end].inbound_types --> [Message{GaussianDistribution}, Message{GaussianDistribution}, Void]
        @fact algo.schedule[end].outbound_type --> GaussianDistribution
        @fact isdefined(algo.schedule[end], :approximation) --> false
    end
end

facts("SumProduct message passing tests") do
    context("SumProduct execute()") do
        context("Should correctly execute a schedule and return the result of the last step") do
            initializeAdditionNode()

            algo = SumProduct(n(:add_node).i[:out])
            prepare!(algo)
            msg = execute(algo)

            @fact msg.payload --> GaussianDistribution(m=0.0, V=2.0)
        end
    end
end
