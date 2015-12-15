facts("VMP.collectInbounds() tests") do
    context("collectInbounds() should add the proper message/marginal") do
        # Mean field factorized Gaussian node
        initializeGaussianNode()
        algo = VMP.Algorithm()
        @fact VMP.collectInbounds(n(:node).i[:mean]) --> (1, [nothing, vague(GammaDistribution), vague(GaussianDistribution)])
        @fact VMP.collectInbounds(n(:node).i[:precision]) --> (2, [vague(GaussianDistribution), nothing, vague(GaussianDistribution)])
        @fact VMP.collectInbounds(n(:node).i[:out]) --> (3, [vague(GaussianDistribution), vague(GammaDistribution), nothing])

        # Structurally factorized
        initializeGaussianNode()
        algo = VMP.Algorithm(Set{Edge}(Edge[n(:node).i[:out].edge])) # Split off extensions of these groups into separate subgraphs
        @fact VMP.collectInbounds(n(:node).i[:mean]) --> (1, [nothing, Message(GammaDistribution()), vague(GaussianDistribution)])
        @fact VMP.collectInbounds(n(:node).i[:precision]) --> (2, [Message(GaussianDistribution()), nothing, vague(GaussianDistribution)])
        @fact VMP.collectInbounds(n(:node).i[:out]) --> (3, [vague(NormalGammaDistribution), vague(NormalGammaDistribution), nothing])
    end
end
