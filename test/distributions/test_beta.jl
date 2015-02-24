#####################
# Unit tests
#####################

facts("BetaDistribution unit tests") do
    context("BetaDistribution() should initiatize a beta distribution") do
        dist = BetaDistribution(a=2.0, b=0.5)
        @fact dist.a => 2.0
        @fact dist.b => 0.5
    end

    context("vague() should initialize a vague (almost uninformative) beta distribution") do
        dist = vague(BetaDistribution)
        @fact dist.a => tiny()
        @fact dist.b => tiny()
    end
end

facts("Marginal calculations for the beta") do

    FactorGraph()

    context("calculateMarginal!(edge) should give correct result and save the marginal to the edge") do
        edge = Edge(TerminalNode(BetaDistribution(a=1.0, b=2.0)),
                    TerminalNode(BetaDistribution(a=3.0, b=4.0)), BetaDistribution)
        calculateForwardMessage!(edge)
        calculateBackwardMessage!(edge)
        marginal_dist = calculateMarginal!(edge)
        @fact edge.marginal => marginal_dist
        @fact edge.marginal.a => 3.0
        @fact edge.marginal.b => 5.0
    end

    context("calculateMarginal(forward_msg, backward_msg) should give correct result") do
        marginal_dist = calculateMarginal(
                                BetaDistribution(a=1.0, b=2.0),
                                BetaDistribution(a=3.0, b=4.0))
        @fact marginal_dist.a => 3.0
        @fact marginal_dist.b => 5.0
    end
end
