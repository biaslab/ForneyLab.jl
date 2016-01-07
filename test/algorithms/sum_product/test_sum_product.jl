facts("collectInbounds() tests") do
    context("collectInbounds() should collect the required inbound messages in an array") do
        # Standard
        initializeGaussianNode()
        @fact ForneyLab.collectInbounds(n(:node).i[:out], Val{symbol(sumProduct!)}) --> (3, [n(:node).i[:mean].partner.message, n(:node).i[:precision].partner.message, nothing])

        # Composite node
        initializeGainEqualityNode(eye(1), Any[Message(DeltaDistribution(1.0)), Message(DeltaDistribution(2.0)), Message(DeltaDistribution(3.0))])
        @fact ForneyLab.collectInbounds(n(:gec_node).i[:out], Val{symbol(sumProduct!)}) --> (3, [n(:gec_node).i[:in1].partner.message, n(:gec_node).i[:in2].partner.message, nothing])
    end
end

facts("SumProduct message passing tests") do
    context("SumProduct execute()") do
        context("Should correctly execute a schedule and return the result of the last step") do
            initializeAdditionNode(Any[Message(GaussianDistribution()), Message(GaussianDistribution()), Message(GaussianDistribution())])

            algo = SumProduct(n(:add_node).i[:out])
            prepare!(algo)
            msg = execute(algo.schedule)

            @fact msg.payload --> GaussianDistribution(m=2.0, V=2.0)
        end

        context("Should correctly execute a schedule and return the result of the last step with preset messages") do
            initializeLoopyGraph(A=[2.0], B=[0.5], noise_m=1.0, noise_V=0.1)
            setMessage!(n(:add).i[:in1], Message(GaussianDistribution(m=2.0, V=0.5)))
            setMessage!(n(:add).i[:out], Message(GaussianDistribution()))
            
            algo = SumProduct(n(:add).i[:in2])
            prepare!(algo)
            msg = execute(algo.schedule)
            dist = ensureParameters!(msg.payload, (:m, :V))

            @fact dist --> n(:add).i[:in2].message.payload
            @fact isApproxEqual(dist.m, [2.0]) --> true
            @fact isApproxEqual(dist.V, reshape([1.5], 1, 1)) --> true
        end

        context("Should handle post-processing of messages (sample)") do
            initializeAdditionNode(Any[Message(GaussianDistribution()), Message(GaussianDistribution()), Message(GaussianDistribution())])

            algo = SumProduct(n(:add_node).i[:out], post_processing_functions=Dict{Interface, Function}(n(:add_node).i[:out] => sample))
            prepare!(algo)
            msg = execute(algo.schedule)

            @fact typeof(msg.payload) --> DeltaDistribution{Float64}
        end
    end
end
