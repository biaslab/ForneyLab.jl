#####################
# Unit tests
#####################

facts("InverseGammaDistribution unit tests") do
    context("InverseGammaDistribution() should initiatize an inverse gamma distribution") do
        dist = InverseGammaDistribution(a=2.0, b=0.5)
        @fact dist.a => 2.0
        @fact dist.b => 0.5
    end

    context("uninformative() should initialize an uninformative inverse gamma distribution") do
        dist = uninformative(InverseGammaDistribution)
        @fact dist.a => -0.999
        @fact dist.b => 0.001
    end
end

facts("Marginal calculations for the inverse gamma") do

    FactorGraph()

    context("calculateMarginal!(edge) should give correct result and save the marginal to the edge") do
        edge = Edge(TerminalNode(InverseGammaDistribution(a=1.0, b=2.0)),
                    TerminalNode(InverseGammaDistribution(a=3.0, b=4.0)), InverseGammaDistribution)
        calculateForwardMessage!(edge)
        calculateBackwardMessage!(edge)
        marginal_dist = calculateMarginal!(edge)
        @fact edge.marginal => marginal_dist
        @fact edge.marginal.a => 5.0
        @fact edge.marginal.b => 6.0
    end

    context("calculateMarginal(forward_msg, backward_msg) should give correct result") do
        marginal_dist = calculateMarginal(
                                InverseGammaDistribution(a=1.0, b=2.0),
                                InverseGammaDistribution(a=3.0, b=4.0))
        @fact marginal_dist.a => 5.0
        @fact marginal_dist.b => 6.0
    end
end