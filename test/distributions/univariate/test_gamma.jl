#####################
# Unit tests
#####################

facts("Gamma unit tests") do
    context("Gamma() should initiatize a gamma distribution") do
        dist = Gamma(a=2.0, b=0.5)
        @fact dist.a --> 2.0
        @fact dist.b --> 0.5
        @fact mean(dist) --> 4.0
        @fact var(dist) --> 8.0
        @fact isnan(mean(Gamma(a=0.0, b=0.5))) --> true
        @fact isnan(mean(Gamma(a=1.0, b=0.0))) --> true
    end

    context("vague() should initialize a vague (almost uninformative) gamma distribution") do
        dist = vague(Gamma)
        @fact dist.a --> tiny
        @fact dist.b --> tiny
    end

    context("calculateMarginal!(edge) should give correct result and save the marginal to the edge") do
        initializePairOfTerminalNodes(Gamma(a=1.0, b=2.0), Gamma(a=3.0, b=4.0))
        edge = n(:t1).i[:out].edge
        n(:t1).i[:out].message = Message(Gamma(a=1.0, b=2.0))
        n(:t2).i[:out].message = Message(Gamma(a=3.0, b=4.0))
        marginal_dist = calculateMarginal!(edge)
        @fact edge.marginal --> marginal_dist
        @fact edge.marginal.a --> 3.0
        @fact edge.marginal.b --> 6.0
    end

    context("Marginal calculation for the combination of a Gamma and Delta") do
        initializePairOfTerminalNodes(Gamma(), Delta(3.0))
        edge = n(:t1).i[:out].edge
        calculateMarginal!(edge)
        @fact edge.marginal --> Delta(3.0)
    end

    context("prod! for Gamma") do
        @fact Gamma(a=1.0, b=2.0) * Gamma(a=3.0, b=4.0) --> Gamma(a=3.0, b=6.0)
        @fact Delta(0.5) * Gamma(a=3.0, b=4.0) --> Delta(0.5)
        @fact_throws DomainError Gamma(a=3.0, b=4.0) * Delta(-1.1)
        @fact typeof(ForneyLab.prod!(Delta(0.5), Gamma(a=3.0, b=4.0), Gamma())) --> Gamma
        @fact mean(ForneyLab.prod!(Gamma(a=3.0, b=4.0), Delta(0.5), Gamma())) --> roughly(0.5)
        @fact var(ForneyLab.prod!(Gamma(a=3.0, b=4.0), Delta(0.5), Gamma())) --> less_than(1e-6)
    end
end
