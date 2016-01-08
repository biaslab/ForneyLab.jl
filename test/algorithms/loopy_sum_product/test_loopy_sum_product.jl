facts("LoopySumProduct message passing tests") do
    context("LoopySumProduct execute()") do
        context("Should correctly execute a schedule and return the result of the last step") do
            initializeAdditionNode([GaussianDistribution(), GaussianDistribution(), GaussianDistribution()])

            algo = LoopySumProduct(n(:add_node).i[:out])
            prepare!(algo)
            msg = execute(algo.schedule)

            @fact msg.payload --> GaussianDistribution(m=2.0, V=2.0)
        end

        context("Should correctly execute a schedule and return the result of the last step with preset messages") do
            initializeLoopyGraph(A=[2.0], B=[0.5], noise_m=1.0, noise_V=0.1)
            
            algo = LoopySumProduct( n(:add).i[:in2],
                                    breaker_messages = Dict(n(:add).i[:in1] => Message(GaussianDistribution(m=2.0, V=0.5)), 
                                                            n(:add).i[:out] => Message(GaussianDistribution())))
            prepare!(algo)
            msg = execute(algo.schedule)
            dist = ensureParameters!(msg.payload, (:m, :V))

            @fact is(dist, n(:add).i[:in2].message.payload) --> true
            @fact isApproxEqual(dist.m, [2.0]) --> true
            @fact isApproxEqual(dist.V, [1.5].') --> true
        end

        context("Should handle post-processing of messages (sample)") do
            initializeAdditionNode([GaussianDistribution(), GaussianDistribution(), GaussianDistribution()])

            algo = LoopySumProduct(n(:add_node).i[:out], post_processing_functions=Dict{Interface, Function}(n(:add_node).i[:out] => sample))
            prepare!(algo)
            msg = execute(algo.schedule)

            @fact typeof(msg.payload) --> DeltaDistribution{Float64}
        end
    end
end