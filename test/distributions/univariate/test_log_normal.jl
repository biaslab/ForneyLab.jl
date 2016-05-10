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

    context("prod!() should yield correct result") do
        @fact DeltaDistribution(0.5) * LogNormalDistribution(m=2.0, s=0.5) --> DeltaDistribution(0.5)
        @fact_throws DomainError LogNormalDistribution(m=2.0, s=0.5) * DeltaDistribution(-1.1)
        @fact typeof(ForneyLab.prod!(DeltaDistribution(0.5), LogNormalDistribution(m=2.0, s=0.5), LogNormalDistribution())) --> LogNormalDistribution
        @fact mean(ForneyLab.prod!(LogNormalDistribution(m=2.0, s=0.5), DeltaDistribution(0.5), LogNormalDistribution())) --> roughly(0.5)
        @fact var(ForneyLab.prod!(LogNormalDistribution(m=2.0, s=0.5), DeltaDistribution(0.5), LogNormalDistribution())) --> less_than(1e-6)
    end
end