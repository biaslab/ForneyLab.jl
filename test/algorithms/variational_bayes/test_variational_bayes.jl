facts("collectInbounds() tests") do
    context("collectInbounds() should add the proper message/marginal") do
        # Mean field factorized Gaussian node
        initializeGaussianNode()
        algo = VariationalBayes()
        @fact ForneyLab.collectInbounds(n(:node).i[:mean], Val{symbol(vmp!)}) --> (1, [nothing, vague(GammaDistribution), vague(GaussianDistribution)])
        @fact ForneyLab.collectInbounds(n(:node).i[:precision], Val{symbol(vmp!)}) --> (2, [vague(GaussianDistribution), nothing, vague(GaussianDistribution)])
        @fact ForneyLab.collectInbounds(n(:node).i[:out], Val{symbol(vmp!)}) --> (3, [vague(GaussianDistribution), vague(GammaDistribution), nothing])

        # Structurally factorized
        initializeGaussianNode()
        algo = VariationalBayes(Set{Edge}(Edge[n(:node).i[:out].edge])) # Split off extensions of these groups into separate subgraphs
        @fact ForneyLab.collectInbounds(n(:node).i[:mean], Val{symbol(vmp!)}) --> (1, [nothing, Message(GammaDistribution()), vague(GaussianDistribution)])
        @fact ForneyLab.collectInbounds(n(:node).i[:precision], Val{symbol(vmp!)}) --> (2, [Message(GaussianDistribution()), nothing, vague(GaussianDistribution)])
        @fact ForneyLab.collectInbounds(n(:node).i[:out], Val{symbol(vmp!)}) --> (3, [vague(NormalGammaDistribution), vague(NormalGammaDistribution), nothing])
    end
end