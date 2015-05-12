facts("Generic Algorithm constructor integration tests") do
    context("Algorithm(::Schedule) should yield Algorithm that executes the passed Schedule") do
        (t_node, gain_node_1, gain_node_2) = initializeChainOfNodes()
        t_node.value = GaussianDistribution(m=5.0, V=2.0)
        schedule = SumProduct.generateSchedule(gain_node_2.out)
        schedule[end].post_processing = mean
        algo = Algorithm(schedule)
        execute(algo)
        @fact gain_node_2.out.message.payload => DeltaDistribution([(mean(t_node.value)*gain_node_1.A*gain_node_2.A)[1]])
    end
end