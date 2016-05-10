facts("LoopySumProduct message passing tests") do
    context("LoopySumProduct execute()") do
        context("Should correctly execute a schedule and return the result of the last step") do
            initializeAdditionNode()

            algo = LoopySumProduct(n(:add_node).i[:out])
            prepare!(algo)
            execute(algo)

            @fact n(:add_node).i[:out].message.payload --> Gaussian(m=0.0, V=2.0)
        end

        context("Should correctly execute a schedule and return the result of the last step with preset messages") do
            FactorGraph()
            GainNode(gain=1.0, id=:g1)
            GainNode(gain=2.0, id=:g2)
            Edge(n(:g1).i[:out], n(:g2).i[:in])
            Edge(n(:g2).i[:out], n(:g1).i[:in])

            algo = LoopySumProduct( n(:g2).i[:out],
                                    breaker_messages=Dict(n(:g2).i[:out] => Message(Gaussian(m=1.0, V=1.0))),
                                    n_iterations=5)
            prepare!(algo)
            execute(algo)
            @fact n(:g2).i[:out].message.payload --> Gaussian(m=32.0, V=1024.0)
        end

    end
end