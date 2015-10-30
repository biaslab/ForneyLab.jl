#####################
# Unit tests
#####################

facts("InverseGammaDistribution unit tests") do
    context("InverseGammaDistribution() should initiatize an inverse gamma distribution") do
        dist = InverseGammaDistribution(a=2.0, b=0.5)
        @fact dist.a --> 2.0
        @fact dist.b --> 0.5
        @fact mean(dist) --> 0.5
        @fact isnan(var(dist)) --> true
        @fact isnan(mean(InverseGammaDistribution(a=1.0, b=0.5))) --> true
        @fact var(InverseGammaDistribution(a=3.0, b=0.5)) --> 1./16
    end

    context("vague() should initialize a vague (almost uninformative) inverse gamma distribution") do
        dist = vague(InverseGammaDistribution)
        @fact dist.a --> tiny
        @fact dist.b --> tiny
    end

    context("Vague inverse gamma distributions should be equal") do
        @fact vague(InverseGammaDistribution) --> vague(InverseGammaDistribution)
    end
end

facts("Marginal calculations for the inverse gamma") do
    context("calculateMarginal!(edge) should give correct result and save the marginal to the edge") do
        initializePairOfTerminalNodes(InverseGammaDistribution(a=1.0, b=2.0), InverseGammaDistribution(a=3.0, b=4.0))
        edge = n(:t1).i[:out].edge
        n(:t1).i[:out].message = Message(InverseGammaDistribution(a=1.0, b=2.0))
        n(:t2).i[:out].message = Message(InverseGammaDistribution(a=3.0, b=4.0))
        marginal_dist = calculateMarginal!(edge)
        @fact edge.marginal --> marginal_dist
        @fact edge.marginal.a --> 5.0
        @fact edge.marginal.b -->  6.0
    end

    context("Marginal calculation for the combination of a InverseGamma and DeltaDistribution") do
        initializePairOfTerminalNodes(InverseGammaDistribution(), DeltaDistribution(3.0))
        edge = n(:t1).i[:out].edge
        calculateMarginal!(edge)
        @fact edge.marginal --> DeltaDistribution(3.0)
    end

    context("calculateMarginal(forward_msg, backward_msg) should give correct result") do
        marginal_dist = calculateMarginal(
                                InverseGammaDistribution(a=1.0, b=2.0),
                                InverseGammaDistribution(a=3.0, b=4.0))
        @fact marginal_dist.a --> 5.0
        @fact marginal_dist.b --> 6.0
    end
end
