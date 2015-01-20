#####################
# Unit tests
#####################

facts("InverseGammaDistribution unit tests") do
    context("InverseGammaDistribution() should initiatize an inverse gamma distribution") do
        dist = InverseGammaDistribution(a=2.0, b=0.5)
        @fact dist.a => 2.0
        @fact dist.b => 0.5
        @fact mean(dist) => 0.5
        @fact isnan(var(dist)) => true
        @fact isnan(mean(InverseGammaDistribution(a=1.0, b=0.5))) => true
        @fact var(InverseGammaDistribution(a=3.0, b=0.5)) => 1./16
    end

    context("vague() should initialize a vague (almost uninformative) inverse gamma distribution") do
        dist = vague(InverseGammaDistribution)
        @fact dist.a => -1.0+tiny()
        @fact dist.b => tiny()
    end

    context("Vague inverse gamma distributions should be equal") do
        @fact vague(InverseGammaDistribution) => vague(InverseGammaDistribution)
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