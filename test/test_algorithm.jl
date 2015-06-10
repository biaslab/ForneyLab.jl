facts("Generic Algorithm constructor integration tests") do
    context("Algorithm(::Schedule) should yield Algorithm that executes the passed Schedule") do
        initializeChainOfNodes()
        n(:node1).value = GaussianDistribution(m=5.0, V=2.0)
        schedule = SumProduct.generateSchedule(n(:node3).i[:out])
        schedule[end].post_processing = mean
        algo = Algorithm(schedule)
        @fact_throws deepcopy(algo)
        execute(algo)
        @fact n(:node3).i[:out].message.payload => DeltaDistribution([(mean(n(:node1).value)*n(:node2).A*n(:node3).A)[1]])
    end
end