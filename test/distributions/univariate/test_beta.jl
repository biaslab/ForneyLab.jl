#####################
# Unit tests
#####################

facts("BetaDistribution unit tests") do
    context("BetaDistribution() should initiatize a beta distribution") do
        dist = BetaDistribution(a=2.0, b=0.5)
        @fact dist.a --> 2.0
        @fact dist.b --> 0.5
    end

    context("vague() should initialize a vague (almost uninformative) beta distribution") do
        dist = vague(BetaDistribution)
        @fact dist.a --> tiny
        @fact dist.b --> tiny
    end
end

facts("Marginal calculations for the beta") do
    context("calculateMarginal!(edge) should give correct result and save the marginal to the edge") do
        initializePairOfTerminalNodes(BetaDistribution(a=1.0, b=2.0), BetaDistribution(a=3.0, b=4.0))
        edge = n(:t1).i[:out].edge
        n(:t1).i[:out].message = Message(BetaDistribution(a=1.0, b=2.0))
        n(:t2).i[:out].message = Message(BetaDistribution(a=3.0, b=4.0))
        marginal_dist = calculateMarginal!(edge)
        @fact edge.marginal --> marginal_dist
        @fact edge.marginal.a --> 3.0
        @fact edge.marginal.b --> 5.0
    end

    context("prod!() should yield correct result") do
        @fact BetaDistribution(a=1.0, b=2.0) * BetaDistribution(a=3.0, b=4.0) --> BetaDistribution(a=3.0, b=5.0)
        @fact DeltaDistribution(0.5) * BetaDistribution(a=3.0, b=4.0) --> DeltaDistribution(0.5)
        @fact_throws DomainError BetaDistribution(a=3.0, b=4.0) * DeltaDistribution(1.1)
        @fact typeof(ForneyLab.prod!(DeltaDistribution(0.5), BetaDistribution(a=3.0, b=4.0), BetaDistribution())) --> BetaDistribution
        @fact mean(ForneyLab.prod!(BetaDistribution(a=3.0, b=4.0), DeltaDistribution(0.5), BetaDistribution())) --> roughly(0.5)
    end
end
