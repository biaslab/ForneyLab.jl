#####################
# Unit tests
#####################

facts("Beta unit tests") do
    context("Beta() should initiatize a beta distribution") do
        dist = Beta(a=2.0, b=0.5)
        @fact dist.a --> 2.0
        @fact dist.b --> 0.5
    end

    context("vague() should initialize a vague (almost uninformative) beta distribution") do
        dist = vague(Beta)
        @fact dist.a --> tiny
        @fact dist.b --> tiny
    end
end

facts("Marginal calculations for the beta") do
    context("calculateMarginal!(edge) should give correct result and save the marginal to the edge") do
        initializePairOfTerminalNodes(Beta(a=1.0, b=2.0), Beta(a=3.0, b=4.0))
        edge = n(:t1).i[:out].edge
        n(:t1).i[:out].message = Message(Beta(a=1.0, b=2.0))
        n(:t2).i[:out].message = Message(Beta(a=3.0, b=4.0))
        marginal_dist = calculateMarginal!(edge)
        @fact edge.marginal --> marginal_dist
        @fact edge.marginal.a --> 3.0
        @fact edge.marginal.b --> 5.0
    end

    context("prod!() should yield correct result") do
        @fact Beta(a=1.0, b=2.0) * Beta(a=3.0, b=4.0) --> Beta(a=3.0, b=5.0)
        @fact Delta(0.5) * Beta(a=3.0, b=4.0) --> Delta(0.5)
        @fact_throws DomainError Beta(a=3.0, b=4.0) * Delta(1.1)
        @fact typeof(ForneyLab.prod!(Delta(0.5), Beta(a=3.0, b=4.0), Beta())) --> Beta
        @fact mean(ForneyLab.prod!(Beta(a=3.0, b=4.0), Delta(0.5), Beta())) --> roughly(0.5)
    end
end
