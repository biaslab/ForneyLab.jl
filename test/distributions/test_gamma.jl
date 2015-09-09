#####################
# Unit tests
#####################

facts("GammaDistribution unit tests") do
    context("GammaDistribution() should initiatize a gamma distribution") do
        dist = GammaDistribution(a=2.0, b=0.5)
        @fact dist.a => 2.0
        @fact dist.b => 0.5
        @fact mean(dist) => 4.0
        @fact var(dist) => 8.0
        @fact isnan(mean(GammaDistribution(a=0.0, b=0.5))) => true
        @fact isnan(mean(GammaDistribution(a=1.0, b=0.0))) => true
    end

    context("vague() should initialize a vague (almost uninformative) gamma distribution") do
        dist = vague(GammaDistribution)
        @fact dist.a => tiny()
        @fact dist.b => tiny()
    end
end

facts("Marginal calculations for the gamma") do
    context("calculateMarginal!(edge) should give correct result and save the marginal to the edge") do
        initializePairOfTerminalNodes(GammaDistribution(a=1.0, b=2.0), GammaDistribution(a=3.0, b=4.0))
        edge = n(:t1).i[:out].edge
        n(:t1).i[:out].message = Message(GammaDistribution(a=1.0, b=2.0))
        n(:t2).i[:out].message = Message(GammaDistribution(a=3.0, b=4.0))
        marginal_dist = calculateMarginal!(edge)
        @fact edge.marginal => marginal_dist
        @fact edge.marginal.a => 3.0
        @fact edge.marginal.b => 6.0
    end

    context("Marginal calculation for the combination of a Gamma and DeltaDistribution") do
        initializePairOfTerminalNodes(GammaDistribution(), DeltaDistribution(3.0))
        edge = n(:t1).i[:out].edge
        calculateMarginal!(edge)
        @fact edge.marginal => DeltaDistribution(3.0)
    end

    context("calculateMarginal(forward_msg, backward_msg) should give correct result") do
        marginal_dist = calculateMarginal(
                                GammaDistribution(a=1.0, b=2.0),
                                GammaDistribution(a=3.0, b=4.0))
        @fact marginal_dist.a => 3.0
        @fact marginal_dist.b => 6.0
    end
end
