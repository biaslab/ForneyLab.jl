#####################
# Unit tests
#####################

facts("InverseGamma unit tests") do
    context("InverseGamma() should initiatize an inverse gamma distribution") do
        dist = InverseGamma(a=2.0, b=0.5)
        @fact dist.a --> 2.0
        @fact dist.b --> 0.5
        @fact mean(dist) --> 0.5
        @fact isnan(var(dist)) --> true
        @fact isnan(mean(InverseGamma(a=1.0, b=0.5))) --> true
        @fact var(InverseGamma(a=3.0, b=0.5)) --> 1./16
    end

    context("vague() should initialize a vague (almost uninformative) inverse gamma distribution") do
        dist = vague(InverseGamma)
        @fact dist.a --> tiny
        @fact dist.b --> tiny
    end

    context("Vague inverse gamma distributions should be equal") do
        @fact vague(InverseGamma) --> vague(InverseGamma)
    end

    context("calculateMarginal!(edge) should give correct result and save the marginal to the edge") do
        initializePairOfTerminalNodes(InverseGamma(a=1.0, b=2.0), InverseGamma(a=3.0, b=4.0))
        edge = n(:t1).i[:out].edge
        n(:t1).i[:out].message = Message(InverseGamma(a=1.0, b=2.0))
        n(:t2).i[:out].message = Message(InverseGamma(a=3.0, b=4.0))
        marginal_dist = calculateMarginal!(edge)
        @fact edge.marginal --> marginal_dist
        @fact edge.marginal.a --> 5.0
        @fact edge.marginal.b -->  6.0
    end

    context("Marginal calculation for the combination of a InverseGamma and Delta") do
        initializePairOfTerminalNodes(InverseGamma(), Delta(3.0))
        edge = n(:t1).i[:out].edge
        calculateMarginal!(edge)
        @fact edge.marginal --> Delta(3.0)
    end

    context("prod!() should yield correct result") do
        @fact InverseGamma(a=1.0, b=2.0) * InverseGamma(a=3.0, b=4.0) --> InverseGamma(a=5.0, b=6.0)
        @fact Delta(0.5) * InverseGamma(a=3.0, b=4.0) --> Delta(0.5)
        @fact_throws DomainError InverseGamma(a=3.0, b=4.0) * Delta(-1.1)
        @fact typeof(ForneyLab.prod!(Delta(0.5), InverseGamma(a=3.0, b=4.0), InverseGamma())) --> InverseGamma
        @fact mean(ForneyLab.prod!(InverseGamma(a=3.0, b=4.0), Delta(0.5), InverseGamma())) --> roughly(0.5)
        @fact var(ForneyLab.prod!(InverseGamma(a=3.0, b=4.0), Delta(0.5), InverseGamma())) --> less_than(1e-6)
    end
end
