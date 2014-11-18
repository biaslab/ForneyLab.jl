#####################
# Unit tests
#####################

facts("GammaDistribution unit tests") do
    context("GammaDistribution() should initiatize a gamma distribution") do
        dist = GammaDistribution(a=2.0, b=0.5)
        @fact dist.a => 2.0
        @fact dist.b => 0.5
    end

    context("uninformative() should initialize an uninformative gamma distribution") do
        dist = uninformative(GammaDistribution)
        @fact dist.a => 1.0-tiny()
        @fact dist.b => tiny()
    end
end

facts("Marginal calculations for the gamma") do

    FactorGraph()

    context("calculateMarginal!(edge) should give correct result and save the marginal to the edge") do
        edge = Edge(TerminalNode(GammaDistribution(a=1.0, b=2.0)),
                    TerminalNode(GammaDistribution(a=3.0, b=4.0)), GammaDistribution)
        calculateForwardMessage!(edge)
        calculateBackwardMessage!(edge)
        marginal_dist = calculateMarginal!(edge)
        @fact edge.marginal => marginal_dist
        @fact edge.marginal.a => 3.0
        @fact edge.marginal.b => 6.0
    end

    context("calculateMarginal(forward_msg, backward_msg) should give correct result") do
        marginal_dist = calculateMarginal(
                                GammaDistribution(a=1.0, b=2.0),
                                GammaDistribution(a=3.0, b=4.0))
        @fact marginal_dist.a => 3.0
        @fact marginal_dist.b => 6.0
    end
end
