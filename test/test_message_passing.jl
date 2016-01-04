#####################
# Unit tests
#####################

facts("Message passing tests") do
    context("execute()") do
        context("Should correctly execute a schedule and return the result of the last step") do
            initializeAdditionNode(Any[Message(GaussianDistribution()), Message(GaussianDistribution()), Message(GaussianDistribution())])

            algo = SumProduct(n(:add_node).i[:out])
            prepare!(algo)
            msg = execute(algo.schedule)

            @fact msg.payload --> GaussianDistribution(m=2.0, V=2.0)
        end

        context("Should correctly execute a schedule and return the result of the last step with preset messages") do
            # TODO: fix

            # initializeLoopyGraph(A=[2.0], B=[0.5], noise_m=1.0, noise_V=0.1)
            # setMessage!(n(:add).i[:in1], Message(GaussianDistribution(m=2.0, V=0.5)))
            # setMessage!(n(:add).i[:out], Message(GaussianDistribution()))
            
            # algo = SumProduct(n(:add).i[:in2])
            # prepare!(algo)
            # msg = execute(algo.schedule)
            # dist = ensureParameters!(msg.payload, (:m, :V))

            # @fact dist --> n(:add).i[:in2].message.payload
            # @fact isApproxEqual(dist.m, [2.0]) --> true
            # @fact isApproxEqual(dist.V, reshape([1.5], 1, 1)) --> true

            @fact true --> false
        end

        context("Should handle post-processing of messages (sample)") do
            initializeAdditionNode(Any[Message(GaussianDistribution()), Message(GaussianDistribution()), Message(GaussianDistribution())])

            algo = SumProduct(n(:add_node).i[:out])
            setPostProcessing!(algo.schedule, n(:add_node).i[:out], sample)
            prepare!(algo)
            msg = execute(algo.schedule)
            dist = ensureParameters!(msg.payload, (:m, :V))

            @fact typeof(dist) --> DeltaDistribution{Float64}
        end
    end

    context("clearMessage!() and clearMessages!() should clear messages") do
        initializeAdditionNode(Any[Message(GaussianDistribution()), Message(GaussianDistribution()), Message(GaussianDistribution())])
        n(:add_node).i[:in1].message = Message(GaussianDistribution())
        n(:add_node).i[:in2].message = Message(GaussianDistribution())
        n(:add_node).i[:out].message = Message(GaussianDistribution())

        clearMessage!(n(:add_node).i[:in2])
        @fact n(:add_node).i[:in2].message --> nothing
        @fact n(:add_node).i[:in1].message --> Message(GaussianDistribution())
        @fact n(:add_node).i[:out].message --> Message(GaussianDistribution())
        clearMessages!(n(:add_node))
        @fact n(:add_node).i[:in1].message --> nothing
        @fact n(:add_node).i[:out].message --> nothing
    end
end
