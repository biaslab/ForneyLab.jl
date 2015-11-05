#####################
# Unit tests
#####################

facts("LogNormalDistribution unit tests") do
    context("LogNormalDistribution() should initiatize a log-normal distribution") do
        dist = LogNormalDistribution(m=2.0, s=0.5)
        @fact dist.m --> 2.0
        @fact dist.s --> 0.5
        @fact mean(dist) --> exp(dist.m + 0.5*dist.s)
        @fact var(dist) --> (exp(dist.s) - 1)*exp(2*dist.m + dist.s)
        @fact isnan(mean(LogNormalDistribution(m=0.0, s=-1.0))) --> true
    end

    context("vague() should initialize a vague (almost uninformative) log-normal distribution") do
        dist = vague(LogNormalDistribution)
        @fact dist.m --> 0.0
        @fact dist.s --> huge
    end
end

facts("Marginal calculations for the log-normal") do
    context("calculateMarginal!(edge) should give correct result and save the marginal to the edge") do
        initializePairOfTerminalNodes(LogNormalDistribution(m=1.0, s=2.0), LogNormalDistribution(m=3.0, s=4.0))
        edge = n(:t1).i[:out].edge
        n(:t1).i[:out].message = Message(LogNormalDistribution(m=1.0, s=2.0))
        n(:t2).i[:out].message = Message(LogNormalDistribution(m=3.0, s=4.0))
        marginal_dist = calculateMarginal!(edge)
        @fact edge.marginal --> marginal_dist
        @fact edge.marginal.m --> 4.0
        @fact edge.marginal.s --> 6.0
    end

    context("Marginal calculation for the combination of a LogNormal and DeltaDistribution") do
        initializePairOfTerminalNodes(LogNormalDistribution(), DeltaDistribution(3.0))
        edge = n(:t1).i[:out].edge
        calculateMarginal!(edge)
        @fact edge.marginal --> DeltaDistribution(3.0)
    end

    context("calculateMarginal(forward_msg, backward_msg) should give correct result") do
        marginal_dist = calculateMarginal(
                                LogNormalDistribution(m=1.0, s=2.0),
                                LogNormalDistribution(m=3.0, s=4.0))
        @fact marginal_dist.m --> 4.0
        @fact marginal_dist.s --> 6.0
    end
end
