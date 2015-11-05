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