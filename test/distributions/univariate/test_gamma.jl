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
        @fact ForneyLab.isProper(Gamma(a=0.0, b=0.5)) --> false
        @fact ForneyLab.isProper(Gamma(a=1.0, b=0.0)) --> false
        @fact pdf(Gamma(a=2.0, b=2.0), 1.5) --> roughly(0.2987224102071837, atol=1e-8)
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
    end

    context("unsafeLogMean() should return correct result") do
        @fact ForneyLab.unsafeLogMean(Gamma(a=1.0, b=2.0)) --> digamma(1) - log(2)
    end

    context("H() should evaluate the entropy") do
        @fact ForneyLab.H(Gamma(a=0.5, b=0.5)) --> roughly(0.7837571104739337)
    end
end
