facts("SumProduct collect inbound type tests") do
    context("SumProduct should collect the proper inbound types") do
        # Standard
        initializeAdditionNode()
        algo = SumProduct(n(:add_node).i[:out])
        @fact algo.schedule[3].inbound_types --> [Message{GaussianDistribution}, Message{GaussianDistribution}, Void]
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
